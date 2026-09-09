# =============================================================================
# DATASET MAESTRO PARA EDA - 2010-2025
#
# Una fila = una celda del soporte aprobado + un año + una semana.
# Conserva las alternativas POWER para compararlas durante el EDA/modelado.
# =============================================================================

library(terra)

ANIOS <- 2010:2025

detectar_raiz <- function() {
  p <- getwd()
  if (requireNamespace("rstudioapi",quietly=TRUE) && rstudioapi::isAvailable()) {
    a <- tryCatch(rstudioapi::getActiveDocumentContext()$path,error=function(e) "")
    if (nzchar(a)) p <- dirname(a)
  }
  repeat {
    if (dir.exists(file.path(p,"datos_proyecto_1"))) return(normalizePath(p))
    q <- dirname(p); if (identical(p,q)) break; p <- q
  }
  stop("No se encontro la raiz. Abra Proyecto1.Rproj.")
}

RAIZ <- detectar_raiz()
BASE <- file.path(RAIZ,"avance_team_proyecto1","reinicio_desde_cero")
DIR_ALT <- file.path(BASE,"power_reconstruido_2010_2025",
                     "alternativas_vecino_bilineal")
DIR_OUT <- file.path(BASE,"dataset_eda")
dir.create(DIR_OUT,recursive=TRUE,showWarnings=FALSE)

s <- readRDS(file.path(BASE,"soporte_resultados","soporte_maestro.rds"))
mascara <- rast(s$mascara); celdas <- s$celdas
xy <- xyFromCell(mascara,celdas)
alt <- rast(file.path(RAIZ,"datos_proyecto_1","imagenes_semanales",
                      "altitud_valle.tif"))
altitud <- values(alt,mat=FALSE)[celdas]

extraer_largo <- function(r,celdas) as.vector(t(values(r,mat=TRUE)[celdas,,drop=FALSE]))

partes <- vector("list",length(ANIOS))
for (i in seq_along(ANIOS)) {
  anio <- ANIOS[i]
  cat("Construyendo año",anio,"...\n")
  precip <- rast(file.path(RAIZ,"datos_proyecto_1","imagenes_semanales",
                           sprintf("chirps_semanal_valle_%d.tif",anio)))
  tn <- rast(file.path(DIR_ALT,sprintf("temperatura_near_%d.tif",anio)))
  tb <- rast(file.path(DIR_ALT,sprintf("temperatura_bilinear_%d.tif",anio)))
  rn <- rast(file.path(DIR_ALT,sprintf("radiacion_near_%d.tif",anio)))
  rb <- rast(file.path(DIR_ALT,sprintf("radiacion_bilinear_%d.tif",anio)))

  # Orden: para cada celda aparecen sus semanas 1...52 consecutivamente.
  partes[[i]] <- data.frame(
    id_celda=rep(seq_along(celdas),each=52),
    celda_raster=rep(celdas,each=52),
    x=rep(xy[,1],each=52), y=rep(xy[,2],each=52),
    altitud_m=rep(altitud,each=52),
    anio=anio, semana=rep(1:52,times=length(celdas)),
    precipitacion_mm=extraer_largo(precip,celdas),
    temperatura_vecino_c=extraer_largo(tn,celdas),
    temperatura_bilineal_c=extraer_largo(tb,celdas),
    radiacion_vecino_mj_m2_dia=extraer_largo(rn,celdas),
    radiacion_bilineal_mj_m2_dia=extraer_largo(rb,celdas)
  )
}
datos_eda <- do.call(rbind,partes)

# CHIRPS utiliza negativos como códigos inválidos. No hay negativos dentro del
# soporte aprobado según auditoría, pero la regla queda explícita y verificable.
datos_eda$precipitacion_mm[datos_eda$precipitacion_mm < 0] <- NA_real_

esperadas <- length(celdas)*52*length(ANIOS)
stopifnot(nrow(datos_eda)==esperadas)

diccionario <- data.frame(
  variable=names(datos_eda),
  descripcion=c(
    "Identificador estable 1-686 del centro dentro del Valle",
    "Número de celda en la plantilla raster 39x37",
    "Longitud del centro, WGS84", "Latitud del centro, WGS84",
    "Altitud SRTM remuestreada a la plantilla, metros",
    "Año calendario", "Bloque semanal de 7 días desde el 1 de enero",
    "Precipitación CHIRPS semanal acumulada, mm",
    "T2M semanal NASA asignada por vecino, °C",
    "T2M semanal NASA remuestreada bilinealmente, °C",
    "Radiación semanal NASA por vecino, MJ/m²/día",
    "Radiación semanal NASA bilineal, MJ/m²/día"),
  rol=c("clave","trazabilidad","coordenada","coordenada","covariable",
        "tiempo","tiempo","respuesta","alternativa","candidata_principal",
        "alternativa","alternativa")
)

cobertura <- data.frame(
  variable=names(datos_eda),
  n_total=nrow(datos_eda),
  n_validos=sapply(datos_eda,function(x) sum(!is.na(x))),
  n_na=sapply(datos_eda,function(x) sum(is.na(x)))
)
cobertura$pct_validos <- 100*cobertura$n_validos/cobertura$n_total

saveRDS(datos_eda,file.path(DIR_OUT,"dataset_eda_2010_2025.rds"),compress=TRUE)
write.csv(diccionario,file.path(DIR_OUT,"diccionario_variables.csv"),row.names=FALSE)
write.csv(cobertura,file.path(DIR_OUT,"control_cobertura.csv"),row.names=FALSE)

# CSV pequeño de ejemplo para inspección sin cargar las 570 mil filas.
set.seed(2026)
muestra <- datos_eda[sample.int(nrow(datos_eda),min(2000,nrow(datos_eda))),]
write.csv(muestra,file.path(DIR_OUT,"muestra_2000_filas.csv"),row.names=FALSE)

cat("\nDATASET EDA CONSTRUIDO\n")
cat(strrep("=",72),"\n",sep="")
cat("Filas:",format(nrow(datos_eda),big.mark=".",decimal.mark=","),"\n")
cat("Celdas:",length(unique(datos_eda$id_celda)),"\n")
cat("Años:",min(datos_eda$anio),"a",max(datos_eda$anio),"\n")
cat("Semanas por año:",length(unique(datos_eda$semana)),"\n\n")
print(cobertura,row.names=FALSE)
cat("\nPrimeras filas:\n")
print(head(datos_eda),row.names=FALSE)
cat("\nPara leerlo en R:\n")
cat('datos_eda <- readRDS("avance_team_proyecto1/reinicio_desde_cero/dataset_eda/dataset_eda_2010_2025.rds")\n')
cat("Directorio:",DIR_OUT,"\n")
