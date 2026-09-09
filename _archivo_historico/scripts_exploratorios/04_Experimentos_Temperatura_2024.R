# =============================================================================
# EXPERIMENTOS FORENSES DE TEMPERATURA NASA POWER - 2024
#
# Busca explicaciones reproducibles para la diferencia con el archivo entregado:
# horario, inicio de semana y alineacion espacial. No modifica datos originales
# ni convierte el mejor ajuste matematico en una correccion automatica.
# =============================================================================

library(terra)

detectar_raiz <- function() {
  candidatos <- c(getwd(), dirname(getwd()), dirname(dirname(getwd())))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error=function(e) "")
    if (nzchar(p)) candidatos <- c(dirname(p), candidatos)
  }
  for (p in candidatos) {
    p <- normalizePath(p, mustWork=FALSE)
    repeat {
      if (dir.exists(file.path(p, "datos_proyecto_1"))) return(p)
      q <- dirname(p); if (identical(p, q)) break; p <- q
    }
  }
  stop("No se encontro la raiz del proyecto.")
}

RAIZ <- detectar_raiz()
BASE <- file.path(RAIZ, "avance_team_proyecto1", "reinicio_desde_cero")
RAW <- file.path(BASE, "power_regional_2024", "raw_nasa")
SALIDA <- file.path(BASE, "power_regional_2024", "resultados_experimentos_temp")
dir.create(SALIDA, recursive=TRUE, showWarnings=FALSE)

soporte <- readRDS(file.path(BASE, "soporte_resultados", "soporte_maestro.rds"))
mascara <- rast(soporte$mascara); plantilla <- mascara; values(plantilla) <- NA
celdas <- soporte$celdas; valle <- vect(soporte$valle)
entregada <- rast(file.path(RAIZ, "datos_proyecto_1", "imagenes_semanales",
                            "power_semanal_valle_2024.tif"))[[1:52]]
obs <- values(entregada, mat=TRUE)[celdas, , drop=FALSE]

semanal <- function(diario, desfase=0L) {
  rast(lapply(1:52, function(s) {
    inicio <- 1L + (s-1L)*7L + desfase
    fin <- inicio + 6L
    if (inicio < 1L || fin > nlyr(diario)) return(diario[[1]] * NA)
    mean(diario[[inicio:fin]], na.rm=TRUE)
  }))
}

evaluar <- function(candidato, experimento, opcion, semanas=2:52) {
  est <- values(candidato[[semanas]], mat=TRUE)[celdas, , drop=FALSE]
  ob <- obs[, semanas, drop=FALSE]
  ok <- is.finite(ob) & is.finite(est); er <- ob[ok] - est[ok]
  data.frame(experimento, opcion, n=sum(ok), ME=mean(er),
             MAE=mean(abs(er)), RMSE=sqrt(mean(er^2)),
             correlacion=cor(ob[ok], est[ok]), max_abs=max(abs(er)))
}

lst_diario <- rast(file.path(RAW, "T2M_2024_LST.nc"))
utc_diario <- rast(file.path(RAW, "T2M_2024_UTC.nc"))
lst_sem <- semanal(lst_diario)

# E1: horario LST contra UTC, manteniendo exactamente lo demas.
res <- rbind(
  evaluar(mask(resample(lst_sem, plantilla, method="bilinear"), mascara),
          "horario", "LST"),
  evaluar(mask(resample(semanal(utc_diario), plantilla, method="bilinear"), mascara),
          "horario", "UTC")
)

# E2: semanas que empiezan 0, 1 o 2 dias despues del 1 de enero.
for (d in 0:2) {
  cand <- mask(resample(semanal(lst_diario, d), plantilla, method="bilinear"), mascara)
  res <- rbind(res, evaluar(cand, "inicio_semana", paste0("dia_1_mas_", d)))
}

# E3: prueba fina de alineacion. Un desplazamiento que mejore solo el RMSE no
# se adopta sin metadatos: sirve para determinar si el error parece geometrico.
desplazamientos <- seq(-0.10, 0.10, by=0.0125)
rejilla <- expand.grid(dx=desplazamientos, dy=desplazamientos)
rejilla$RMSE <- rejilla$MAE <- rejilla$ME <- NA_real_
for (i in seq_len(nrow(rejilla))) {
  movida <- shift(lst_sem, dx=rejilla$dx[i], dy=rejilla$dy[i])
  cand <- mask(resample(movida, plantilla, method="bilinear"), mascara)
  z <- evaluar(cand, "alineacion", sprintf("dx_%+.4f_dy_%+.4f",
                                           rejilla$dx[i], rejilla$dy[i]))
  rejilla[i, c("RMSE", "MAE", "ME")] <- z[1, c("RMSE", "MAE", "ME")]
}
mejor <- rejilla[which.min(rejilla$RMSE), ]
base <- rejilla[rejilla$dx == 0 & rejilla$dy == 0, ]
res <- rbind(res, data.frame(
  experimento="alineacion", opcion=sprintf("mejor_dx_%+.4f_dy_%+.4f",
                                            mejor$dx, mejor$dy),
  n=NA, ME=mejor$ME, MAE=mejor$MAE, RMSE=mejor$RMSE,
  correlacion=NA, max_abs=NA))

# E4: una calibracion lineal se calcula solo como diagnostico de versionado.
base_r <- mask(resample(lst_sem, plantilla, method="bilinear"), mascara)
est <- values(base_r, mat=TRUE)[celdas, 2:52, drop=FALSE]
ob <- obs[, 2:52, drop=FALSE]; ok <- is.finite(ob) & is.finite(est)
modelo_cal <- lm(ob[ok] ~ est[ok])
pred <- predict(modelo_cal); er <- ob[ok] - pred
ecuacion <- sprintf("entregada = %.6f + %.6f * nasa",
                    coef(modelo_cal)[1], coef(modelo_cal)[2])
res <- rbind(res, data.frame(experimento="calibracion_diagnostica",
  opcion=ecuacion, n=sum(ok), ME=mean(er), MAE=mean(abs(er)),
  RMSE=sqrt(mean(er^2)), correlacion=cor(ob[ok], pred), max_abs=max(abs(er))))

write.csv(res, file.path(SALIDA, "comparacion_experimentos.csv"), row.names=FALSE)
write.csv(rejilla, file.path(SALIDA, "rejilla_desplazamientos.csv"), row.names=FALSE)

cat("\nEXPERIMENTOS TEMPERATURA 2024 — COMPARACION S02-S52\n")
cat(strrep("=", 76), "\n", sep="")
print(res, row.names=FALSE, digits=6)
cat("\nAlineacion base: RMSE", round(base$RMSE, 6), "\n")
cat("Mejor desplazamiento exploratorio: dx", mejor$dx, "dy", mejor$dy,
    "RMSE", round(mejor$RMSE, 6), "\n")
cat("Un desplazamiento optimizado NO se aprobara sin metadatos de procedencia.\n")

dibujar <- function() {
  par(mfrow=c(2,2), mar=c(4.5,4.5,3,1))
  z <- res[res$experimento == "horario",]
  barplot(z$RMSE, names.arg=z$opcion, col=c("steelblue","tan"),
          ylab="RMSE (°C)", main="Horario")
  z <- res[res$experimento == "inicio_semana",]
  barplot(z$RMSE, names.arg=c("1 ene","2 ene","3 ene"), col="steelblue",
          ylab="RMSE (°C)", main="Inicio de bloques semanales")
  m <- matrix(rejilla$RMSE, nrow=length(desplazamientos), byrow=FALSE)
  image(desplazamientos, desplazamientos, m, xlab="Desplazamiento lon (°)",
        ylab="Desplazamiento lat (°)", main="RMSE por alineación",
        col=hcl.colors(30, "YlOrRd", rev=TRUE))
  points(mejor$dx, mejor$dy, pch=4, lwd=3)
  plot(est[ok], ob[ok], pch=16, cex=.25, col=rgb(0,0,0,.2),
       xlab="NASA actual bilineal (°C)", ylab="Entregada (°C)",
       main="Calibración solo diagnóstica")
  abline(0,1,lty=2,col="red"); abline(modelo_cal,col="blue",lwd=2)
}
dibujar()
pdf(file.path(SALIDA, "experimentos_temperatura_2024.pdf"), 11, 8.5)
dibujar(); dev.off()
cat("Resultados:", SALIDA, "\n")
