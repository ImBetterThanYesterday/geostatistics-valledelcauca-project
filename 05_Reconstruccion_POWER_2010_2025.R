# =============================================================================
# RECONSTRUCCION NASA POWER 2010-2025
#
# Metodos aprobados provisionalmente por evidencia 2024:
# - T2M: promedio semanal y remuestreo bilineal.
# - ALLSKY_SFC_SW_DWN: promedio semanal y vecino mas cercano.
# Conserva NetCDF crudos, no sobrescribe datos entregados y valida cada año.
# =============================================================================

library(terra)

ANIOS <- 2010:2025
BBOX <- c(lon_min=-78, lon_max=-75.5, lat_min=2.5, lat_max=5.5)

detectar_raiz <- function() {
  p <- getwd()
  if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
    a <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error=function(e) "")
    if (nzchar(a)) p <- dirname(a)
  }
  repeat {
    if (dir.exists(file.path(p, "datos_proyecto_1"))) return(normalizePath(p))
    q <- dirname(p); if (identical(p,q)) break; p <- q
  }
  stop("No se encontro la raiz. Abra Proyecto1.Rproj.")
}

RAIZ <- detectar_raiz()
BASE <- file.path(RAIZ, "avance_team_proyecto1", "reinicio_desde_cero")
DIR_PROCESO <- file.path(BASE, "power_reconstruido_2010_2025")
DIR_RAW <- file.path(DIR_PROCESO, "raw_nasa")
DIR_TIF <- file.path(DIR_PROCESO, "geotiff_semanal")
DIR_RES <- file.path(DIR_PROCESO, "resultados_validacion")
for (d in c(DIR_RAW, DIR_TIF, DIR_RES)) dir.create(d, recursive=TRUE,
                                                   showWarnings=FALSE)

soporte <- readRDS(file.path(BASE, "soporte_resultados", "soporte_maestro.rds"))
mascara <- rast(soporte$mascara); plantilla <- mascara; values(plantilla) <- NA
celdas <- soporte$celdas

url_regional <- function(parametro, anio) {
  paste0("https://power.larc.nasa.gov/api/temporal/daily/regional?",
         "parameters=", parametro, "&community=AG",
         "&longitude-min=", BBOX["lon_min"], "&longitude-max=", BBOX["lon_max"],
         "&latitude-min=", BBOX["lat_min"], "&latitude-max=", BBOX["lat_max"],
         "&start=", anio, "0101&end=", anio, "1231",
         "&format=NETCDF&time-standard=LST")
}

descargar_seguro <- function(url, destino, intentos=4L) {
  if (file.exists(destino) && file.info(destino)$size > 1000) {
    cat("  Crudo existente:", basename(destino), "\n"); return(invisible(TRUE))
  }
  for (i in seq_len(intentos)) {
    cat("  Descarga", i, ":", basename(destino), "\n")
    bien <- tryCatch({
      download.file(url, destino, mode="wb", quiet=TRUE)
      file.exists(destino) && file.info(destino)$size > 1000
    }, error=function(e) FALSE)
    if (bien) return(invisible(TRUE))
    if (i < intentos) Sys.sleep(2*i)
  }
  stop("No se pudo descargar despues de ", intentos, " intentos: ", url)
}

agregar_52 <- function(diario) {
  if (nlyr(diario) < 364) stop("Serie diaria incompleta: ", sources(diario)[1])
  r <- rast(lapply(1:52, function(s) {
    dias <- ((s-1)*7+1):(s*7)
    mean(diario[[dias]], na.rm=TRUE)
  }))
  names(r) <- sprintf("semana_%02d", 1:52)
  r
}

metrica <- function(obs, est, variable, anio, semanas=1:52, periodo="S01_S52") {
  o <- values(obs[[semanas]], mat=TRUE)[celdas,,drop=FALSE]
  e <- values(est[[semanas]], mat=TRUE)[celdas,,drop=FALSE]
  ok <- is.finite(o) & is.finite(e); er <- o[ok]-e[ok]
  data.frame(anio, variable, periodo, n=sum(ok), ME=mean(er),
             MAE=mean(abs(er)), RMSE=sqrt(mean(er^2)),
             correlacion=if(sum(ok)>2) cor(o[ok],e[ok]) else NA_real_,
             max_abs=max(abs(er)), coincidencias_1e_4=sum(abs(er)<=1e-4),
             pct_coincidencia=100*mean(abs(er)<=1e-4))
}

validacion <- data.frame()
manifiesto <- data.frame()

for (anio in ANIOS) {
  cat("\n", strrep("-", 65), "\nAÑO ", anio, "\n", sep="")
  f_t <- file.path(DIR_RAW, sprintf("T2M_%d_LST.nc", anio))
  f_r <- file.path(DIR_RAW, sprintf("RAD_%d_LST.nc", anio))
  descargar_seguro(url_regional("T2M", anio), f_t)
  descargar_seguro(url_regional("ALLSKY_SFC_SW_DWN", anio), f_r)

  d_t <- rast(f_t); d_r <- rast(f_r)
  if (nlyr(d_t) != ifelse(anio %% 4 == 0, 366, 365))
    stop("Numero inesperado de dias T2M en ", anio)
  if (nlyr(d_r) != nlyr(d_t)) stop("T2M y radiacion difieren en dias: ", anio)

  temp <- mask(resample(agregar_52(d_t), plantilla, method="bilinear"), mascara)
  rad <- mask(resample(agregar_52(d_r), plantilla, method="near"), mascara)
  names(temp) <- sprintf("temp_semana_%02d", 1:52)
  names(rad) <- sprintf("radiacion_semana_%02d", 1:52)
  combinado <- c(temp, rad)
  f_out <- file.path(DIR_TIF, sprintf("power_reconstruido_valle_%d.tif", anio))
  writeRaster(combinado, f_out, overwrite=TRUE, datatype="FLT4S")

  entregado <- rast(file.path(RAIZ, "datos_proyecto_1", "imagenes_semanales",
                              sprintf("power_semanal_valle_%d.tif", anio)))
  validacion <- rbind(validacion,
    metrica(entregado[[1:52]], temp, "temperatura", anio),
    metrica(entregado[[1:52]], temp, "temperatura", anio, 2:52, "S02_S52"),
    metrica(entregado[[53:104]], rad, "radiacion", anio),
    metrica(entregado[[53:104]], rad, "radiacion", anio, 2:52, "S02_S52"))
  manifiesto <- rbind(manifiesto, data.frame(
    anio, dias=nlyr(d_t), inicio=as.character(time(d_t)[1]),
    fin=as.character(tail(time(d_t),1)), temp_crudo=basename(f_t),
    rad_crudo=basename(f_r), salida=basename(f_out), celdas_soporte=length(celdas)))
  print(tail(validacion, 4), row.names=FALSE, digits=5)
}

write.csv(validacion, file.path(DIR_RES, "validacion_por_anio.csv"), row.names=FALSE)
write.csv(manifiesto, file.path(DIR_RES, "manifiesto_procedencia.csv"), row.names=FALSE)

dibujar <- function() {
  par(mfrow=c(2,2), mar=c(4.5,4.5,3,1))
  for (v in c("temperatura", "radiacion")) {
    z <- validacion[validacion$variable==v & validacion$periodo=="S02_S52",]
    plot(z$anio, z$RMSE, type="b", pch=19, xlab="Año", ylab="RMSE",
         main=paste("Validación", v, "S02-S52"))
    abline(h=0, lty=2, col="grey40")
  }
  z <- validacion[validacion$variable=="radiacion" & validacion$periodo=="S02_S52",]
  barplot(z$pct_coincidencia, names.arg=z$anio, las=2, cex.names=.7,
          ylim=c(0,100), ylab="Coincidencia exacta (%)",
          main="Radiación: réplica por año S02-S52")
  z <- validacion[validacion$variable=="temperatura" & validacion$periodo=="S02_S52",]
  plot(z$anio, z$correlacion, type="b", pch=19, ylim=c(min(z$correlacion)-.005,1),
       xlab="Año", ylab="Correlación", main="Temperatura: correlación anual")
}
dibujar()
pdf(file.path(DIR_RES, "resumen_validacion_2010_2025.pdf"), 11, 8.5)
dibujar(); dev.off()

cat("\nRECONSTRUCCION COMPLETADA\n")
cat("NetCDF preservados:", DIR_RAW, "\n")
cat("GeoTIFF candidatos:", DIR_TIF, "\n")
cat("Validacion:", DIR_RES, "\n")
