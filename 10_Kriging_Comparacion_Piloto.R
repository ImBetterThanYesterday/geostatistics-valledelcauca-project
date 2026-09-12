# =============================================================================
# 10_Kriging_Comparacion_Piloto.R
#
# Compara kriging ordinario y universal para el corte 2014-banda 38.
# Se usa la tendencia cuadrática ambiental y el semivariograma esférico
# seleccionado provisionalmente. Validación LOO con 686 celdas.
# =============================================================================
suppressPackageStartupMessages({library(sp); library(gstat); library(terra)})

buscar_base <- function() {
  cc <- unique(normalizePath(c(getwd(),file.path(getwd(),".."),
    file.path(getwd(),"../.."),file.path(getwd(),"../../..")),mustWork=FALSE))
  ok <- vapply(cc,function(p) file.exists(file.path(p,
    "avance_team_proyecto1","reinicio_desde_cero","tendencia_flexible_resultados",
    "residuales_cuadratica_ambiental.csv")),logical(1))
  if(!any(ok)) stop("No se encontró el residual cuadrático.")
  cc[which(ok)[1]]
}
BASE <- buscar_base()
CUT <- readRDS(file.path(BASE,"avance_team_proyecto1","reinicio_desde_cero",
                         "cortes_espaciales_resultados","tres_cortes_espaciales.rds"))$intermedia
CUT <- CUT[order(CUT$id_celda),]
lat0 <- mean(CUT$y)
CUT$x_km <- (CUT$x-mean(CUT$x))*111.32*cos(lat0*pi/180)
CUT$y_km <- (CUT$y-lat0)*111.32
CUT$x2_km <- CUT$x_km^2
CUT$y2_km <- CUT$y_km^2
CUT$xy_km <- CUT$x_km*CUT$y_km
CUT$.y <- sqrt(CUT$precipitacion_mm)

resid <- read.csv(file.path(BASE,"avance_team_proyecto1","reinicio_desde_cero",
                             "tendencia_flexible_resultados",
                             "residuales_cuadratica_ambiental.csv"))
vr <- read.csv(file.path(BASE,"avance_team_proyecto1","reinicio_desde_cero",
                         "semivariograma_cuadratico_resultados",
                         "comparacion_modelos_variograma_cuadratico.csv"))
p <- vr[vr$modelo=="esferico",]
if(nrow(p)!=1 || !is.finite(p$range_km)) stop("No hay parámetros esféricos válidos.")
modelo_vgm <- vgm(psill=p$psill, model="Sph", range=p$range_km, nugget=p$nugget)

# Fórmulas: ordinario estima una media constante; universal usa la tendencia.
f_ordinario <- .y ~ 1
f_universal <- .y ~ x_km + y_km + x2_km + y2_km + xy_km +
  altitud_m + temperatura_bilineal_c + radiacion_bilineal_mj_m2_dia

sp_all <- SpatialPointsDataFrame(coords=CUT[,c("x_km","y_km")],
                                  data=CUT[,c(".y","x_km","y_km","x2_km","y2_km",
                                    "xy_km","altitud_m","temperatura_bilineal_c",
                                    "radiacion_bilineal_mj_m2_dia")])
predicciones <- list()
for(tipo in c(ordinario="ordinario", universal="universal")) {
  f <- if(tipo=="ordinario") f_ordinario else f_universal
  pred <- numeric(nrow(CUT)); varp <- numeric(nrow(CUT))
  for(i in seq_len(nrow(CUT))) {
    train <- sp_all[-i,]
    nuevo <- sp_all[i,]
    k <- tryCatch(krige(f, locations=train, newdata=nuevo, model=modelo_vgm,
                        nmax=100), error=function(e) NULL)
    if(is.null(k)) { pred[i] <- NA; varp[i] <- NA } else {
      pred[i] <- pmax(k$var1.pred,0)^2
      varp[i] <- k$var1.var
    }
  }
  e <- CUT$precipitacion_mm-pred
  predicciones[[tipo]] <- data.frame(id_celda=CUT$id_celda,x=CUT$x,y=CUT$y,
    observado_mm=CUT$precipitacion_mm,predicho_mm=pred,error_mm=e,
    varianza_kriging=varp,tipo_kriging=tipo)
}
metricas <- do.call(rbind,lapply(predicciones,function(d) {
  ok <- is.finite(d$predicho_mm); e <- d$error_mm[ok]
  data.frame(tipo_kriging=unique(d$tipo_kriging), n_validos=sum(ok),
    MAE=mean(abs(e)), RMSE=sqrt(mean(e^2)), sesgo=mean(-e),
    R2_predictivo=1-sum(e^2)/sum((d$observado_mm[ok]-mean(d$observado_mm[ok]))^2))
}))
OUT <- file.path(BASE,"avance_team_proyecto1","reinicio_desde_cero",
                 "kriging_piloto_resultados")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
write.csv(metricas,file.path(OUT,"comparacion_kriging_LOO.csv"),row.names=FALSE)
for(nm in names(predicciones))
  write.csv(predicciones[[nm]],file.path(OUT,paste0("predicciones_",nm,"_LOO.csv")),row.names=FALSE)
print(metricas,row.names=FALSE)

# Predicción sobre la superficie completa para ambos métodos.
grid <- sp_all
superficies <- list()
for(tipo in c(ordinario="ordinario", universal="universal")) {
  f <- if(tipo=="ordinario") f_ordinario else f_universal
  k <- krige(f, locations=sp_all, newdata=grid, model=modelo_vgm, nmax=100)
  superficies[[tipo]] <- data.frame(x=CUT$x,y=CUT$y,
    predicho_raiz=k$var1.pred, varianza=k$var1.var)
  write.csv(superficies[[tipo]],file.path(OUT,paste0("superficie_",tipo,".csv")),
            row.names=FALSE)
}

pdf(file.path(OUT,"comparacion_kriging_ordinario_universal.pdf"),width=11,height=8.5)
par(mfrow=c(2,2),mar=c(4.5,4.5,3.3,1.5))
for(tipo in names(predicciones)) {
  d<-predicciones[[tipo]]
  plot(d$observado_mm,d$predicho_mm,pch=19,cex=.6,col=adjustcolor("steelblue",.6),
    xlab="Observado (mm)",ylab="Predicho LOO (mm)",main=paste("LOO:",tipo))
  abline(0,1,lty=2,lwd=2)
}
for(tipo in names(predicciones)) {
  d<-predicciones[[tipo]]
  r<-rast(d[,c("x","y","error_mm")],type="xyz",crs="EPSG:4326")
  plot(r,col=hcl.colors(40,"RdYlBu",rev=TRUE),
    main=paste("Error LOO:",tipo),xlab="Longitud",ylab="Latitud")
}
dev.off()
cat("Resultados guardados en:",OUT,"\n")
