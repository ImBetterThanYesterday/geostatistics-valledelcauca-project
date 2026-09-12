#######################################################################
# PUNTO 0 — CIERRE DE DECISIONES ANTES DEL MODELO DE MEDIA
# No ajusta modelos, Moran, variogramas ni kriging.
#######################################################################

library(terra)

buscar_raiz <- function() {
  p <- normalizePath(getwd(), mustWork=TRUE)
  repeat {
    if (dir.exists(file.path(p,"datos_proyecto_1")) &&
        dir.exists(file.path(p,"avance_team_proyecto1"))) return(p)
    padre <- dirname(p); if (identical(p,padre)) break; p <- padre
  }
  stop("No se encontro Proyecto1. Abra la carpeta raiz en RStudio.")
}

RAIZ <- buscar_raiz()
BASE <- file.path(RAIZ,"avance_team_proyecto1","reinicio_desde_cero")
DIR_RAW <- file.path(BASE,"power_reconstruido_2010_2025","raw_nasa")
DIR_DATOS <- file.path(RAIZ,"datos_proyecto_1")
DIR_OUT <- file.path(BASE,"punto_0_resultados")
dir.create(DIR_OUT,recursive=TRUE,showWarnings=FALSE)

cat("========== PUNTO 0: SOPORTE TEMPORAL Y ESCENARIOS ==========\n")

# 1. Integridad diaria NASA POWER.
anios <- 2010:2025
es_bisiesto <- function(a) (a %% 4 == 0 & a %% 100 != 0) | a %% 400 == 0
integridad <- do.call(rbind,lapply(anios,function(a) {
  ft <- file.path(DIR_RAW,sprintf("T2M_%d_LST.nc",a))
  fr <- file.path(DIR_RAW,sprintf("RAD_%d_LST.nc",a))
  if (!file.exists(ft) || !file.exists(fr)) stop("Faltan NetCDF de ",a)
  rt <- rast(ft); rr <- rast(fr); tt <- time(rt)
  esperados <- ifelse(es_bisiesto(a),366L,365L)
  data.frame(anio=a,dias_temperatura=nlyr(rt),dias_radiacion=nlyr(rr),
             fecha_inicio=as.character(min(tt)),fecha_fin=as.character(max(tt)),
             dias_esperados=esperados,serie_completa=nlyr(rt)==esperados && nlyr(rr)==esperados,
             dias_fuera_de_52_bloques=esperados-364L)
}))
print(integridad,row.names=FALSE)

# 2. Evidencia de la regla semanal con CHIRPS 2024.
fs_chirps <- file.path(DIR_DATOS,"cache_raster",sprintf("chirps_diario_2024_%02d.tif",1:12))
if (!all(file.exists(fs_chirps))) stop("Faltan archivos diarios CHIRPS 2024")
diario_2024 <- rast(fs_chirps)
semanal_2024 <- rast(file.path(DIR_DATOS,"imagenes_semanales","chirps_semanal_valle_2024.tif"))
soporte <- readRDS(file.path(BASE,"soporte_resultados","soporte_maestro.rds"))
celdas_soporte <- soporte$celdas
banda_52_7dias <- sum(diario_2024[[358:364]],na.rm=TRUE)
banda_52_9dias <- sum(diario_2024[[358:366]],na.rm=TRUE)

comparar <- function(obs,est) {
  if(!compareGeom(obs,est,stopOnError=FALSE)) est <- resample(est,obs,method="bilinear")
  x <- values(obs,mat=FALSE)[celdas_soporte]
  y <- values(est,mat=FALSE)[celdas_soporte]
  ok <- is.finite(x)&is.finite(y)
  c(n=sum(ok),MAE=mean(abs(x[ok]-y[ok])),correlacion=cor(x[ok],y[ok]))
}
m7 <- comparar(semanal_2024[[52]],banda_52_7dias)
m9 <- comparar(semanal_2024[[52]],banda_52_9dias)
prueba_regla <- data.frame(
  regla=c("Dias 358-364 (7 dias)","Dias 358-366 (9 dias)"),
  n=c(m7["n"],m9["n"]), MAE=c(m7["MAE"],m9["MAE"]),
  correlacion=c(m7["correlacion"],m9["correlacion"]), row.names=NULL
)
cat("\nComparacion de la banda 52 de CHIRPS 2024:\n")
print(prueba_regla,row.names=FALSE)

# 3. Dataset y escenarios que pasan a modelamiento.
dataset <- readRDS(file.path(BASE,"dataset_eda","dataset_final_7_variables.rds"))
control_dataset <- data.frame(
  criterio=c("filas","celdas","anios","bandas_por_anio","nulos_totales"),
  valor=c(nrow(dataset),length(unique(dataset$id_celda)),length(unique(dataset$anio)),
          length(unique(dataset$semana)),sum(is.na(dataset)))
)

decisiones <- data.frame(
  tema=c("Unidad temporal","Dias 365/366","Cobertura 2025","Temperatura POWER",
         "Radiacion POWER","Unidad de analisis","Respuesta","Validacion"),
  decision=c("52 bloques de 7 dias desde el 1 de enero",
             "Excluir: quedan fuera del esquema de 364 dias entregado",
             "Completa: 2025-01-01 a 2025-12-31",
             "Comparar vecino y bilinear; bilinear es candidata principal",
             "Comparar vecino y bilinear sin ganador previo",
             "Una celda x anio x banda semanal",
             "Comparar original, raiz cuadrada y log1p",
             "Pendiente aprobar particiones temporal, espacial y espacio-temporal"),
  estado=c(rep("APROBADA",7),"PENDIENTE")
)
cat("\nDecisiones del punto 0:\n"); print(decisiones,row.names=FALSE)

write.csv(integridad,file.path(DIR_OUT,"integridad_diaria_2010_2025.csv"),row.names=FALSE)
write.csv(prueba_regla,file.path(DIR_OUT,"prueba_regla_banda_52_chirps_2024.csv"),row.names=FALSE)
write.csv(control_dataset,file.path(DIR_OUT,"control_dataset_modelamiento.csv"),row.names=FALSE)
write.csv(decisiones,file.path(DIR_OUT,"decisiones_punto_0.csv"),row.names=FALSE)

writeLines(c(
  "PUNTO 0 — CONCLUSION",
  "Los NetCDF POWER 2010-2025 estan completos por dia.",
  "Los productos entregados usan 52 bloques fijos de 7 dias (364 dias).",
  "El remanente anual no se imputara ni se agregara a la banda 52.",
  "Vecino y bilinear se compararan como escenarios con las mismas particiones.",
  "No se calculara el variograma definitivo antes de seleccionar el modelo de media.",
  "Pendiente: aprobar el diseno exacto de validacion."
),file.path(DIR_OUT,"conclusion_punto_0.txt"))

cat("\nPunto 0 terminado. Resultados:",DIR_OUT,"\n")
