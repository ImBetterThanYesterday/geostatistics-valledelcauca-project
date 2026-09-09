###################################################################
# EDA por variable - dataset final
# Visualización en RStudio y resultados reproducibles en archivos
###################################################################

library(terra)

buscar_raiz <- function() {
  candidatos <- unique(normalizePath(c(getwd(), file.path(getwd(), ".."),
                                        file.path(getwd(), "../..")), mustWork=FALSE))
  ok <- vapply(candidatos, function(p) dir.exists(file.path(p, "datos_proyecto_1")) &&
                 dir.exists(file.path(p, "avance_team_proyecto1")), logical(1))
  if (!any(ok)) stop("No se encontro la raiz. Abra Proyecto1 en RStudio.")
  candidatos[which(ok)[1]]
}
RAIZ_PROYECTO <- buscar_raiz()
DIR_SCRIPT <- file.path(RAIZ_PROYECTO, "avance_team_proyecto1", "reinicio_desde_cero")
DIR_RESULTADOS <- file.path(DIR_SCRIPT, "eda_final_resultados")
dir.create(DIR_RESULTADOS, recursive=TRUE, showWarnings=FALSE)
dataset <- readRDS(file.path(DIR_SCRIPT, "dataset_eda", "dataset_final_7_variables.rds"))
str(dataset)

variables <- c(
  "precipitacion_mm", "temperatura_bilineal_c",
  "radiacion_bilineal_mj_m2_dia", "altitud_m",
  "clim_precipitacion_mm", "clim_temperatura_bilineal_c",
  "clim_radiacion_bilineal_mj_m2_dia"
)

nombres <- c(
  "Precipitación", "Temperatura", "Radiación solar", "Altitud",
  "Climatología de precipitación", "Climatología de temperatura",
  "Climatología de radiación"
)

unidades <- c(
  "mm por periodo", "°C", "MJ/m²/día", "m",
  "mm por periodo", "°C", "MJ/m²/día"
)

nombres_matriz <- c(
  "Precipitación", "Temperatura", "Radiación", "Altitud",
  "Clim. precip.", "Clim. temp.", "Clim. rad."
)

colores <- c(
  "steelblue", "firebrick", "orange", "darkgreen",
  "steelblue4", "firebrick4", "darkorange3"
)

################################################################
# 0) Cobertura temporal aproximada
################################################################

# Esta aproximación no verifica las fechas reales de las bandas.
dataset$mes_aprox <- pmin(12, ceiling(dataset$semana * 7 / 30.4))
tabla_meses <- table(dataset$anio, dataset$mes_aprox) > 0

print(tabla_meses)
print(rowSums(tabla_meses))

################################################################
# 1) Estadísticas descriptivas y revisión de calidad
################################################################

# Rangos orientativos: estar fuera señala un valor para revisar.
rango_min <- c(0, -10, 0, -50, 0, -10, 0)
rango_max <- c(700, 45, 45, 6000, 700, 45, 45)

for (i in seq_along(variables)) {
  v <- variables[i]
  # Evitar inflar n, atipicos y negativos por repeticion: altitud se resume en
  # 686 celdas; climatologias en 686 x 52; variables dinamicas en todas las filas.
  base_resumen <- if (v == "altitud_m") {
    dataset[!duplicated(dataset$id_celda), ]
  } else if (grepl("^clim_", v)) {
    dataset[!duplicated(dataset[, c("id_celda", "semana")]), ]
  } else dataset
  x <- base_resumen[[v]]
  
  cat("\n====================================================\n")
  cat(nombres[i], "|", unidades[i], "\n")
  cat("====================================================\n")
  print(summary(x))
  
  cat("Desviación estándar:", sd(x, na.rm = TRUE), "\n")
  cat("Nulos:", sum(is.na(x)), "\n")
  cat("Negativos:", sum(x < 0, na.rm = TRUE), "\n")
  cat("Ceros:", sum(x == 0, na.rm = TRUE),
      sprintf("(%.2f%%)\n", 100 * mean(x == 0, na.rm = TRUE)))
  
  fuera_rango <- sum(
    x < rango_min[i] | x > rango_max[i],
    na.rm = TRUE
  )
  cat("Fuera del rango orientativo [", rango_min[i], ",",
      rango_max[i], "]:", fuera_rango, "\n")
  
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  ric <- q[2] - q[1]
  atipicos <- x < q[1] - 1.5 * ric | x > q[2] + 1.5 * ric
  
  cat("Atípicos de Tukey:", sum(atipicos, na.rm = TRUE),
      sprintf("(%.2f%%)\n", 100 * mean(atipicos, na.rm = TRUE)))
  
  if (v != "precipitacion_mm") {
    # Para variables reducidas se compara contra precipitacion en el mismo
    # nivel de agregacion, no contra las 570752 filas repetidas.
    if (v == "altitud_m") {
      p_comp <- aggregate(precipitacion_mm ~ id_celda, dataset, mean)$precipitacion_mm
    } else if (grepl("^clim_", v)) {
      p_comp <- base_resumen$clim_precipitacion_mm
    } else p_comp <- dataset$precipitacion_mm
    r <- cor(x, p_comp, use = "complete.obs")
    cat("Pearson con precipitación en nivel comparable:", round(r, 3), "\n")
  }
}

################################################################
# 2) Matriz de correlación en consola
################################################################

matriz_cor <- cor(dataset[, variables], use = "pairwise.complete.obs")
print(round(matriz_cor, 3))
write.csv(matriz_cor, file.path(DIR_RESULTADOS, "correlacion_pearson_pooled.csv"))

################################################################
# 2.1) Anomalias respecto a la climatologia propia 2010-2025
################################################################

# Una anomalia positiva indica un valor superior al promedio historico de esa
# misma celda y banda semanal. No es una derivada matematica.
dataset$anom_precipitacion_mm <- dataset$precipitacion_mm - dataset$clim_precipitacion_mm
dataset$anom_temperatura_c <- dataset$temperatura_bilineal_c - dataset$clim_temperatura_bilineal_c
dataset$anom_radiacion_mj_m2_dia <- dataset$radiacion_bilineal_mj_m2_dia -
  dataset$clim_radiacion_bilineal_mj_m2_dia

variables_anom <- c("anom_precipitacion_mm", "anom_temperatura_c",
                    "anom_radiacion_mj_m2_dia")

cat("\n========== ANOMALIAS CLIMATICAS ==========" , "\n")
print(summary(dataset[, variables_anom]))
matriz_cor_anom <- cor(dataset[, variables_anom], use="complete.obs")
cat("\nCorrelacion entre anomalias (variacion respecto a lo esperado):\n")
print(round(matriz_cor_anom, 3))

################################################################
# 2.2) Correlaciones en tres escalas
################################################################

vars_base <- c("precipitacion_mm", "temperatura_bilineal_c",
               "radiacion_bilineal_mj_m2_dia", "altitud_m")

# Escala espacial: cada celda aporta una observacion (promedio 2010-2025).
datos_espaciales <- aggregate(dataset[, vars_base],
                              by=list(id_celda=dataset$id_celda), FUN=mean)
cor_espacial <- cor(datos_espaciales[, vars_base], use="complete.obs")

# Escala temporal regional: una observacion por anio y banda semanal.
datos_temporales <- aggregate(dataset[, vars_base[1:3]],
                              by=list(anio=dataset$anio, semana=dataset$semana), FUN=mean)
cor_temporal <- cor(datos_temporales[, vars_base[1:3]], use="complete.obs")

cat("\n========== CORRELACION ESPACIAL (686 promedios por celda) ==========\n")
print(round(cor_espacial, 3))
cat("\n========== CORRELACION TEMPORAL REGIONAL (832 bandas) ==========\n")
print(round(cor_temporal, 3))

cor_tres_escalas <- rbind(
  data.frame(escala="Pooled espacio-tiempo", predictor=vars_base[-1],
             correlacion=matriz_cor["precipitacion_mm", vars_base[-1]]),
  data.frame(escala="Espacial: promedio por celda", predictor=vars_base[-1],
             correlacion=cor_espacial["precipitacion_mm", vars_base[-1]]),
  data.frame(escala="Temporal: promedio regional", predictor=vars_base[2:3],
             correlacion=cor_temporal["precipitacion_mm", vars_base[2:3]]),
  data.frame(escala="Anomalias", predictor=variables_anom[2:3],
             correlacion=matriz_cor_anom["anom_precipitacion_mm", variables_anom[2:3]])
)
print(cor_tres_escalas, row.names=FALSE)

################################################################
# 2.3) Transformaciones candidatas de precipitacion
################################################################

asimetria <- function(x) {
  x <- x[is.finite(x)]; m <- mean(x); s <- sd(x)
  mean((x-m)^3) / s^3
}

transformaciones <- data.frame(
  transformacion=c("Original", "Raiz cuadrada", "log1p"),
  media=c(mean(dataset$precipitacion_mm), mean(sqrt(dataset$precipitacion_mm)),
          mean(log1p(dataset$precipitacion_mm))),
  sd=c(sd(dataset$precipitacion_mm), sd(sqrt(dataset$precipitacion_mm)),
       sd(log1p(dataset$precipitacion_mm))),
  asimetria=c(asimetria(dataset$precipitacion_mm),
              asimetria(sqrt(dataset$precipitacion_mm)),
              asimetria(log1p(dataset$precipitacion_mm)))
)
cat("\n========== TRANSFORMACIONES DE PRECIPITACION ==========\n")
print(transformaciones, row.names=FALSE)
cat("Menor asimetria absoluta no selecciona por si sola el modelo; la decision final\n")
cat("se tomara comparando residuales y validacion predictiva.\n")

################################################################
# 2.4) Colinealidad preliminar entre covariables
################################################################

predictores_vif <- c("altitud_m", "temperatura_bilineal_c",
                     "radiacion_bilineal_mj_m2_dia", "clim_precipitacion_mm",
                     "clim_temperatura_bilineal_c", "clim_radiacion_bilineal_mj_m2_dia")

# VIF calculado con regresiones auxiliares. Se usa una muestra reproducible para
# evitar ajustar seis regresiones sobre filas espacio-temporales repetidas.
set.seed(2026)
idx_vif <- sample.int(nrow(dataset), min(100000L, nrow(dataset)))
datos_vif <- dataset[idx_vif, predictores_vif]
vif <- sapply(predictores_vif, function(v) {
  otros <- setdiff(predictores_vif, v)
  ajuste <- lm(reformulate(otros, response=v), data=datos_vif)
  1 / (1-summary(ajuste)$r.squared)
})
tabla_vif <- data.frame(variable=names(vif), VIF=as.numeric(vif),
                        diagnostico=ifelse(vif>=10,"Alto",ifelse(vif>=5,"Revisar","Aceptable")))
cat("\n========== VIF PRELIMINAR ==========\n")
print(tabla_vif, row.names=FALSE)
cat("El VIF es diagnostico, no una regla automatica para eliminar variables.\n")

################################################################
# 3) Panel de precipitación
################################################################

fila_max <- dataset[which.max(dataset$precipitacion_mm), ]
anio_max <- fila_max$anio
semana_max <- fila_max$semana

sub_max <- dataset[
  dataset$anio == anio_max & dataset$semana == semana_max,
  c("x", "y", "precipitacion_mm")
]

r_max <- rast(sub_max, type = "xyz", crs = "EPSG:4326")

promedio_celda <- aggregate(
  precipitacion_mm ~ id_celda + x + y,
  data = dataset, FUN = mean, na.rm = TRUE
)

r_precip <- rast(
  promedio_celda[, c("x", "y", "precipitacion_mm")],
  type = "xyz", crs = "EPSG:4326"
)

par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2),
    oma = c(0, 0, 2, 0), cex.main = 0.95)

hist(dataset$precipitacion_mm, breaks = 40, col = colores[1],
     border = "white", main = "Histograma",
     xlab = unidades[1], ylab = "Frecuencia")

boxplot(dataset$precipitacion_mm, col = colores[1],
        main = "Boxplot", ylab = unidades[1])

# Márgenes explícitos para evitar el recorte del título.
plot(r_precip,
     main = "Distribución espacial\npromedio",
     mar = c(4.5, 4.5, 4.5, 5),
     cex.main = 0.95,
     xlab = "Longitud", ylab = "Latitud",
     col = hcl.colors(50, "Blues", rev = TRUE))

# Se usa "banda" porque semana es la etiqueta del dataset,
# cuya correspondencia con la semana calendario requiere revisión.
plot(r_max,
     main = paste0("Evento máximo\nBanda ", semana_max, " | ", anio_max),
     mar = c(4.5, 4.5, 4.5, 5),
     cex.main = 0.95,
     xlab = "Longitud", ylab = "Latitud",
     col = hcl.colors(50, "Blues", rev = TRUE))

points(fila_max$x, fila_max$y, pch = 4, col = "red3", lwd = 2)

legend("bottomleft",
       legend = paste(round(fila_max$precipitacion_mm, 2), "mm"),
       pch = 4, col = "red3", bty = "n", cex = 0.8)

mtext("Precipitación | mm por periodo", outer = TRUE,
      side = 3, line = 0.4, font = 2)

cat("\nRegistro del máximo de precipitación:\n")
print(
  fila_max[, c("id_celda", "anio", "semana", "x", "y",
               "precipitacion_mm")],
  row.names = FALSE
)

################################################################
# 4) Panel individual de los seis predictores
################################################################

for (i in 2:length(variables)) {
  v <- variables[i]
  x <- dataset[[v]]
  
  # Relación resumida mediante 20 intervalos del predictor.
  ok <- is.finite(x) & is.finite(dataset$precipitacion_mm)
  bins <- cut(x[ok], breaks = 20, include.lowest = TRUE)
  
  promedio_por_bin <- tapply(dataset$precipitacion_mm[ok], bins, mean)
  centro_bin <- tapply(x[ok], bins, mean)
  bins_validos <- is.finite(centro_bin) & is.finite(promedio_por_bin)
  r <- cor(x[ok], dataset$precipitacion_mm[ok])
  
  datos_mapa <- data.frame(
    id_celda = dataset$id_celda,
    x = dataset$x, y = dataset$y, valor = x
  )
  
  promedio_celda <- aggregate(
    valor ~ id_celda + x + y,
    data = datos_mapa, FUN = mean, na.rm = TRUE
  )
  
  r_v <- rast(
    promedio_celda[, c("x", "y", "valor")],
    type = "xyz", crs = "EPSG:4326"
  )
  
  # Título del mapa según el tipo de variable.
  titulo_mapa <- if (v == "altitud_m") {
    "Altitud del terreno"
  } else if (grepl("^clim_", v)) {
    "Promedio climatológico\npor celda"
  } else {
    "Distribución espacial\npromedio"
  }
  
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2),
      oma = c(0, 0, 2, 0), cex.main = 0.95)
  
  hist(x, breaks = 40, col = colores[i], border = "white",
       main = "Histograma", xlab = unidades[i], ylab = "Frecuencia")
  
  boxplot(x, col = colores[i],
          main = "Boxplot", ylab = unidades[i])
  
  plot(centro_bin[bins_validos], promedio_por_bin[bins_validos],
       type = "o", pch = 16, cex = 0.6, col = colores[i], lwd = 2,
       xlab = unidades[i], ylab = "Precipitación media (mm/periodo)",
       main = paste0("Medias por intervalos\nPearson r = ", round(r, 3)))
  grid()
  
  # Márgenes explícitos para los mapas de todos los predictores.
  plot(r_v,
       main = titulo_mapa,
       mar = c(4.5, 4.5, 4.5, 5),
       cex.main = 0.95,
       xlab = "Longitud", ylab = "Latitud",
       col = hcl.colors(50, "Blues", rev = TRUE))
  
  mtext(paste(nombres[i], "|", unidades[i]), outer = TRUE,
        side = 3, line = 0.4, font = 2)
}

################################################################
# 5) Matriz de correlación como mapa de calor
################################################################

par(mfrow = c(1, 1), mar = c(8, 9, 4, 2),
    oma = c(0, 0, 0, 0), cex.main = 1)

n <- ncol(matriz_cor)
paleta <- colorRampPalette(c("#B2182B", "white", "#2166AC"))(201)

# Invertir filas para mostrar la matriz como una tabla.
image(1:n, 1:n, t(matriz_cor[n:1, ]), zlim = c(-1, 1),
      col = paleta, axes = FALSE, xlab = "", ylab = "",
      main = "Matriz de correlación\nPearson",
      xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5),
      xaxs = "i", yaxs = "i")

axis(1, at = 1:n, labels = nombres_matriz, las = 2,
     tick = FALSE, cex.axis = 0.85)

axis(2, at = 1:n, labels = rev(nombres_matriz), las = 2,
     tick = FALSE, cex.axis = 0.85)

abline(v = seq(0.5, n + 0.5), h = seq(0.5, n + 0.5),
       col = "grey85", lwd = 0.7)

for (fila in 1:n) {
  for (columna in 1:n) {
    valor <- matriz_cor[fila, columna]
    
    color_texto <- if (is.finite(valor) && abs(valor) >= 0.6) {
      "white"
    } else {
      "black"
    }
    
    text(columna, n - fila + 1,
         labels = if (is.finite(valor)) sprintf("%.2f", valor) else "NA",
         col = color_texto, cex = 0.9, font = 2)
  }
}

box()

################################################################
# Restablecer configuración gráfica
################################################################

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1),
    oma = c(0, 0, 0, 0), cex.main = 1.2)

cat("\nEDA terminado. Consulta los paneles anteriores con las flechas de Plots.\n")

# Tablas auditables. Tukey describe valores extremos; no implica eliminarlos.
resumen_eda <- do.call(rbind, lapply(seq_along(variables), function(i) {
  v <- variables[i]
  d <- if (v == "altitud_m") dataset[!duplicated(dataset$id_celda), ] else
    if (grepl("^clim_", v)) dataset[!duplicated(dataset[,c("id_celda","semana")]), ] else dataset
  x <- d[[v]]; q <- quantile(x, c(.25,.75), na.rm=TRUE); ric <- diff(q)
  data.frame(variable=v, unidad=unidades[i], n=length(x), n_na=sum(is.na(x)),
             media=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE), minimo=min(x,na.rm=TRUE),
             p25=q[1], mediana=median(x,na.rm=TRUE), p75=q[2], maximo=max(x,na.rm=TRUE),
             n_negativos=sum(x<0,na.rm=TRUE), n_ceros=sum(x==0,na.rm=TRUE),
             n_extremos_tukey=sum(x<q[1]-1.5*ric | x>q[2]+1.5*ric,na.rm=TRUE),
             nivel_analisis=if(v=="altitud_m") "686 celdas" else if(grepl("^clim_",v)) "686 celdas x 52 semanas" else "espacio-tiempo completo")
}))
write.csv(resumen_eda, file.path(DIR_RESULTADOS, "estadisticas_descriptivas.csv"), row.names=FALSE)
write.csv(as.data.frame.matrix(tabla_meses), file.path(DIR_RESULTADOS, "cobertura_meses_aproximados.csv"))
write.csv(fila_max[,c("id_celda","anio","semana","x","y","precipitacion_mm")],
          file.path(DIR_RESULTADOS, "evento_maximo_precipitacion.csv"), row.names=FALSE)
write.csv(matriz_cor_anom, file.path(DIR_RESULTADOS, "correlacion_anomalias.csv"))
write.csv(cor_espacial, file.path(DIR_RESULTADOS, "correlacion_espacial_promedios.csv"))
write.csv(cor_temporal, file.path(DIR_RESULTADOS, "correlacion_temporal_regional.csv"))
write.csv(cor_tres_escalas, file.path(DIR_RESULTADOS, "correlaciones_tres_escalas.csv"), row.names=FALSE)
write.csv(transformaciones, file.path(DIR_RESULTADOS, "transformaciones_precipitacion.csv"), row.names=FALSE)
write.csv(tabla_vif, file.path(DIR_RESULTADOS, "vif_preliminar.csv"), row.names=FALSE)

# Paneles nuevos visibles en RStudio.
dibujar_transformaciones <- function() {
  par(mfrow=c(1,3), mar=c(4.5,4.2,3,1))
  hist(dataset$precipitacion_mm, breaks=40, col="steelblue", border="white",
       main=paste0("Original\nasimetria = ",round(transformaciones$asimetria[1],2)), xlab="mm")
  hist(sqrt(dataset$precipitacion_mm), breaks=40, col="darkgreen", border="white",
       main=paste0("Raiz cuadrada\nasimetria = ",round(transformaciones$asimetria[2],2)), xlab="sqrt(mm)")
  hist(log1p(dataset$precipitacion_mm), breaks=40, col="darkorange", border="white",
       main=paste0("log1p\nasimetria = ",round(transformaciones$asimetria[3],2)), xlab="log(1 + mm)")
}

dibujar_anomalias <- function() {
  par(mfrow=c(1,3), mar=c(4.5,4.2,3,1))
  hist(dataset$anom_precipitacion_mm, breaks=40, col="steelblue", border="white",
       main="Anomalia precipitacion", xlab="mm frente a climatologia")
  hist(dataset$anom_temperatura_c, breaks=40, col="firebrick", border="white",
       main="Anomalia temperatura", xlab="grados C frente a climatologia")
  hist(dataset$anom_radiacion_mj_m2_dia, breaks=40, col="orange", border="white",
       main="Anomalia radiacion", xlab="MJ/m2/dia frente a climatologia")
}

dibujar_transformaciones()
dibujar_anomalias()
par(mfrow=c(1,1))

# Tres visuales de síntesis en PNG; los paneles completos siguen visibles en RStudio.
png(file.path(DIR_RESULTADOS, "precipitacion_resumen.png"), width=1600, height=1200, res=160)
par(mfrow=c(2,2), mar=c(4.5,4.5,3,3))
hist(dataset$precipitacion_mm, breaks=40, col=colores[1], border="white", main="Precipitacion", xlab=unidades[1])
boxplot(dataset$precipitacion_mm, col=colores[1], main="Boxplot", ylab=unidades[1])
plot(r_precip, main="Promedio espacial", col=hcl.colors(50,"Blues",rev=TRUE))
plot(r_max, main=paste("Maximo: banda",semana_max,"-",anio_max), col=hcl.colors(50,"Blues",rev=TRUE)); points(fila_max$x,fila_max$y,pch=4,col="red",lwd=2)
dev.off()

png(file.path(DIR_RESULTADOS, "matriz_correlacion.png"), width=1500, height=1300, res=170)
par(mar=c(9,10,4,2)); image(1:n,1:n,t(matriz_cor[n:1,]),zlim=c(-1,1),col=paleta,axes=FALSE,xlab="",ylab="",main="Correlacion Pearson (pooled)")
axis(1,at=1:n,labels=nombres_matriz,las=2,tick=FALSE); axis(2,at=1:n,labels=rev(nombres_matriz),las=2,tick=FALSE)
for(fila in 1:n) for(columna in 1:n) text(columna,n-fila+1,sprintf("%.2f",matriz_cor[fila,columna]),col=if(abs(matriz_cor[fila,columna])>=.6) "white" else "black",font=2)
box(); dev.off()

png(file.path(DIR_RESULTADOS, "transformaciones_precipitacion.png"), width=1800, height=650, res=160)
dibujar_transformaciones(); dev.off()
png(file.path(DIR_RESULTADOS, "anomalias_climaticas.png"), width=1800, height=650, res=160)
dibujar_anomalias(); dev.off()

writeLines(c("EDA FINAL - RESUMEN", paste("Fecha:",Sys.time()),
             paste("Dataset:",nrow(dataset),"filas x",ncol(dataset),"columnas"),
             "Soporte: 686 celdas x 16 anios x 52 bandas", "Nulos: 0 en las variables analizadas",
             "Nota: las correlaciones pooled mezclan variacion espacial y temporal; no prueban causalidad.",
             "Se agregaron correlaciones espacial, temporal y de anomalias.",
             paste("Transformacion con menor asimetria absoluta:", transformaciones$transformacion[which.min(abs(transformaciones$asimetria))]),
             paste("Predictores con VIF >= 5:", paste(tabla_vif$variable[tabla_vif$VIF>=5], collapse=", ")),
             "Nota: transformacion y variables definitivas se decidiran con residuales y validacion.",
             "Nota: los extremos de Tukey no se eliminan automaticamente, especialmente en precipitacion."),
           file.path(DIR_RESULTADOS,"resumen_eda.txt"))
cat("Resultados guardados en:", DIR_RESULTADOS, "\n")
