# =============================================================================
# COMPARACION VECINO VS BILINEAL PARA TEMPERATURA Y RADIACION (2010-2025)
#
# Conserva ambos metodos, los compara exclusivamente donde el profesor entrego
# valores validos y crea candidatos separados. No reemplaza archivos originales.
# =============================================================================

library(terra)

ANIOS <- 2010:2025
METODOS <- c("near", "bilinear")

detectar_raiz <- function() {
  p <- getwd()
  if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
    a <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error=function(e) "")
    if (nzchar(a)) p <- dirname(a)
  }
  repeat {
    if (dir.exists(file.path(p,"datos_proyecto_1"))) return(normalizePath(p))
    q <- dirname(p); if (identical(p,q)) break; p <- q
  }
  stop("No se encontro la raiz del proyecto.")
}

RAIZ <- detectar_raiz()
BASE <- file.path(RAIZ,"avance_team_proyecto1","reinicio_desde_cero")
DIR_POWER <- file.path(BASE,"power_reconstruido_2010_2025")
DIR_RAW <- file.path(DIR_POWER,"raw_nasa")
DIR_ALT <- file.path(DIR_POWER,"alternativas_vecino_bilineal")
DIR_RES <- file.path(DIR_POWER,"comparacion_metodos")
dir.create(DIR_ALT,recursive=TRUE,showWarnings=FALSE)
dir.create(DIR_RES,recursive=TRUE,showWarnings=FALSE)

s <- readRDS(file.path(BASE,"soporte_resultados","soporte_maestro.rds"))
mascara <- rast(s$mascara); plantilla <- mascara; values(plantilla) <- NA
celdas <- s$celdas

agregar_52 <- function(x) {
  r <- rast(lapply(1:52,function(sem) {
    d <- ((sem-1)*7+1):(sem*7)
    mean(x[[d]],na.rm=TRUE)
  }))
  names(r) <- sprintf("semana_%02d",1:52)
  r
}

calcular <- function(obs,est,anio,variable,metodo,semanas,periodo) {
  o <- values(obs[[semanas]],mat=TRUE)[celdas,,drop=FALSE]
  e <- values(est[[semanas]],mat=TRUE)[celdas,,drop=FALSE]
  ok <- is.finite(o)&is.finite(e); er <- o[ok]-e[ok]
  data.frame(anio,variable,metodo,periodo,n=sum(ok),ME=mean(er),
             MAE=mean(abs(er)),RMSE=sqrt(mean(er^2)),
             correlacion=cor(o[ok],e[ok]),max_abs=max(abs(er)),
             pct_exacta=100*mean(abs(er)<=1e-4))
}

resultados <- data.frame()
for (anio in ANIOS) {
  cat("Procesando",anio,"...\n")
  t_sem <- agregar_52(rast(file.path(DIR_RAW,sprintf("T2M_%d_LST.nc",anio))))
  r_sem <- agregar_52(rast(file.path(DIR_RAW,sprintf("RAD_%d_LST.nc",anio))))
  prof <- rast(file.path(RAIZ,"datos_proyecto_1","imagenes_semanales",
                         sprintf("power_semanal_valle_%d.tif",anio)))
  for (metodo in METODOS) {
    temp <- mask(resample(t_sem,plantilla,method=metodo),mascara)
    rad <- mask(resample(r_sem,plantilla,method=metodo),mascara)
    names(temp) <- sprintf("temp_%s_semana_%02d",metodo,1:52)
    names(rad) <- sprintf("radiacion_%s_semana_%02d",metodo,1:52)
    writeRaster(temp,file.path(DIR_ALT,sprintf("temperatura_%s_%d.tif",metodo,anio)),
                overwrite=TRUE,datatype="FLT4S")
    writeRaster(rad,file.path(DIR_ALT,sprintf("radiacion_%s_%d.tif",metodo,anio)),
                overwrite=TRUE,datatype="FLT4S")
    resultados <- rbind(resultados,
      calcular(prof[[1:52]],temp,anio,"temperatura",metodo,1:52,"S01_S52"),
      calcular(prof[[1:52]],temp,anio,"temperatura",metodo,2:52,"S02_S52"),
      calcular(prof[[53:104]],rad,anio,"radiacion",metodo,1:52,"S01_S52"),
      calcular(prof[[53:104]],rad,anio,"radiacion",metodo,2:52,"S02_S52"))
  }
}

write.csv(resultados,file.path(DIR_RES,"metricas_por_anio_y_metodo.csv"),row.names=FALSE)

# Resumen global ponderado: se vuelven a acumular errores, evitando promediar
# RMSE anuales como si todos tuvieran distinto peso.
resumen <- aggregate(cbind(MAE,RMSE,correlacion,pct_exacta)~variable+metodo+periodo,
                     resultados,median)
names(resumen)[4:7] <- paste0(names(resumen)[4:7],"_mediana_anual")
write.csv(resumen,file.path(DIR_RES,"resumen_mediana_anual.csv"),row.names=FALSE)

ganadores <- do.call(rbind,lapply(split(resultados,
  interaction(resultados$anio,resultados$variable,resultados$periodo)),function(z) {
    z <- z[which.min(z$RMSE),]
    z[,c("anio","variable","periodo","metodo","RMSE","MAE","correlacion")]
}))
write.csv(ganadores,file.path(DIR_RES,"ganador_por_anio.csv"),row.names=FALSE)

cat("\nRESUMEN: MEDIANA DE METRICAS ANUALES\n")
cat(strrep("=",72),"\n",sep="")
print(resumen,row.names=FALSE,digits=5)
cat("\nMETODO CON MENOR RMSE POR AÑO (S02-S52)\n")
print(table(subset(ganadores,periodo=="S02_S52")[,c("variable","metodo")]))

dibujar <- function() {
  par(mfrow=c(2,2),mar=c(4.5,4.5,3,1))
  for (v in c("temperatura","radiacion")) {
    z <- resultados[resultados$variable==v & resultados$periodo=="S02_S52",]
    rango <- range(z$RMSE)
    plot(NA,xlim=range(ANIOS),ylim=rango,xlab="Año",ylab="RMSE",
         main=paste("RMSE",v,"S02-S52"))
    for (m in METODOS) {
      q <- z[z$metodo==m,]
      lines(q$anio,q$RMSE,type="b",pch=ifelse(m=="near",16,17),
            col=ifelse(m=="near","firebrick","steelblue"))
    }
    legend("topright",c("Vecino","Bilineal"),col=c("firebrick","steelblue"),
           pch=c(16,17),lty=1,bty="n")
  }
  z <- subset(resumen,periodo=="S02_S52")
  for (v in c("temperatura","radiacion")) {
    q <- z[z$variable==v,]
    barplot(q$RMSE_mediana_anual,names.arg=q$metodo,col=c("steelblue","firebrick"),
            ylab="Mediana RMSE anual",main=paste("Comparación global",v))
  }
}
dibujar()
pdf(file.path(DIR_RES,"comparacion_vecino_bilineal_2010_2025.pdf"),11,8.5)
dibujar(); dev.off()
cat("\nAlternativas guardadas en:",DIR_ALT,"\n")
cat("Resultados guardados en:",DIR_RES,"\n")
