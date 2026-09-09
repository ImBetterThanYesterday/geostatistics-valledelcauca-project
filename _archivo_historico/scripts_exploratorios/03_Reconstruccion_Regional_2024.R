# =============================================================================
# PROYECTO 1 - RECONSTRUCCION REGIONAL NASA POWER 2024
#
# Proposito:
#   1. Reconstruir 52 semanas desde los NetCDF diarios originales de NASA.
#   2. Comparar remuestreo vecino, bilineal y cubico contra el GeoTIFF entregado.
#   3. Diagnosticar por que temperatura no coincide exactamente.
#   4. Crear candidatos completos sobre las 686 celdas aprobadas.
#
# Este script NO sobrescribe los datos suministrados y NO corrige el sesgo de
# temperatura. Los productos son candidatos auditables, no datos finales.
# =============================================================================

library(terra)

ANIO <- 2024L
METODOS <- c("near", "bilinear", "cubicspline")

detectar_raiz_proyecto <- function() {
  inicios <- getwd()
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    activo <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                       error = function(e) "")
    if (nzchar(activo)) inicios <- c(inicios, dirname(activo))
  }
  archivo <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (length(archivo) && nzchar(archivo)) inicios <- c(inicios, dirname(archivo))
  for (inicio in unique(inicios)) {
    actual <- normalizePath(inicio, mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(actual, "datos_proyecto_1"))) return(actual)
      padre <- dirname(actual)
      if (identical(padre, actual)) break
      actual <- padre
    }
  }
  stop("No se encontro el proyecto. Abra Proyecto1.Rproj y vuelva a ejecutar.")
}

RAIZ <- detectar_raiz_proyecto()
setwd(RAIZ)
DIR_REINICIO <- file.path(RAIZ, "avance_team_proyecto1", "reinicio_desde_cero")
DIR_REGIONAL <- file.path(DIR_REINICIO, "power_regional_2024")
DIR_RAW <- file.path(DIR_REGIONAL, "raw_nasa")
DIR_SALIDA <- file.path(DIR_REGIONAL, "resultados")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

arch_temp <- file.path(DIR_RAW, "T2M_2024_LST.nc")
arch_rad <- file.path(DIR_RAW, "RAD_2024_LST.nc")
if (!file.exists(arch_temp) || !file.exists(arch_rad)) {
  stop("Faltan los NetCDF diarios en power_regional_2024/raw_nasa.")
}

soporte <- readRDS(file.path(DIR_REINICIO, "soporte_resultados",
                             "soporte_maestro.rds"))
mascara <- rast(soporte$mascara)
valle <- vect(soporte$valle)
celdas_valle <- soporte$celdas
plantilla <- mascara
values(plantilla) <- NA_real_

entregado <- rast(file.path(RAIZ, "datos_proyecto_1", "imagenes_semanales",
                            sprintf("power_semanal_valle_%d.tif", ANIO)))
if (nlyr(entregado) != 104L) stop("POWER entregado debe tener 104 bandas.")
temp_entregada <- entregado[[1:52]]
rad_entregada <- entregado[[53:104]]
names(temp_entregada) <- names(rad_entregada) <- sprintf("semana_%02d", 1:52)

agregar_semanal <- function(r_diario) {
  # Semana del proyecto = 7 dias consecutivos desde el 1 de enero.
  # Se usan dias 1:364; en 2024 quedan fuera 30 y 31 de diciembre.
  if (nlyr(r_diario) < 364L) stop("La serie diaria tiene menos de 364 dias.")
  salidas <- lapply(1:52, function(s) {
    dias <- ((s - 1L) * 7L + 1L):(s * 7L)
    mean(r_diario[[dias]], na.rm = TRUE)
  })
  r <- rast(salidas)
  names(r) <- sprintf("semana_%02d", 1:52)
  r
}

aplicar_soporte <- function(r) {
  salida <- mask(r, mascara)
  names(salida) <- sprintf("semana_%02d", 1:52)
  salida
}

temp_nativa <- agregar_semanal(rast(arch_temp))
rad_nativa <- agregar_semanal(rast(arch_rad))

candidatos_temp <- setNames(lapply(METODOS, function(m) {
  aplicar_soporte(resample(temp_nativa, plantilla, method = m))
}), METODOS)
candidatos_rad <- setNames(lapply(METODOS, function(m) {
  aplicar_soporte(resample(rad_nativa, plantilla, method = m))
}), METODOS)

metricas <- function(observado, estimado, variable, metodo, semanas = 1:52,
                     periodo = "semanas_01_52") {
  o <- values(observado[[semanas]], mat = TRUE)[celdas_valle, , drop = FALSE]
  e <- values(estimado[[semanas]], mat = TRUE)[celdas_valle, , drop = FALSE]
  ok <- is.finite(o) & is.finite(e)
  error <- o[ok] - e[ok]
  data.frame(variable, metodo, periodo, n = length(error),
             ME = mean(error), MAE = mean(abs(error)),
             RMSE = sqrt(mean(error^2)),
             correlacion = cor(o[ok], e[ok]),
             diferencia_maxima = max(abs(error)))
}

tabla_metricas <- do.call(rbind, lapply(METODOS, function(m) rbind(
  metricas(temp_entregada, candidatos_temp[[m]], "temperatura", m),
  metricas(temp_entregada, candidatos_temp[[m]], "temperatura", m, 2:52,
           "semanas_02_52"),
  metricas(rad_entregada, candidatos_rad[[m]], "radiacion", m),
  metricas(rad_entregada, candidatos_rad[[m]], "radiacion", m, 2:52,
           "semanas_02_52")
)))
write.csv(tabla_metricas, file.path(DIR_SALIDA, "metricas_metodos_2024.csv"),
          row.names = FALSE)

# Diagnostico fino del mejor candidato de temperatura: entregado - NASA actual.
o_temp <- values(temp_entregada, mat = TRUE)[celdas_valle, , drop = FALSE]
e_temp <- values(candidatos_temp$bilinear, mat = TRUE)[celdas_valle, , drop = FALSE]
error_temp <- o_temp - e_temp
error_temp[!is.finite(o_temp) | !is.finite(e_temp)] <- NA_real_
sesgo_celda <- rowMeans(error_temp, na.rm = TRUE)
sesgo_celda[!is.finite(sesgo_celda)] <- NA_real_
sesgo_semana <- colMeans(error_temp, na.rm = TRUE)
sd_espacial_semana <- apply(error_temp, 2, sd, na.rm = TRUE)

diagnostico_semanal <- data.frame(
  semana = 1:52, sesgo_medio = sesgo_semana,
  sd_espacial = sd_espacial_semana,
  n = colSums(is.finite(error_temp))
)
diagnostico_celdas <- data.frame(
  celda = celdas_valle,
  x = xyFromCell(plantilla, celdas_valle)[, 1],
  y = xyFromCell(plantilla, celdas_valle)[, 2],
  sesgo_medio = sesgo_celda,
  n = rowSums(is.finite(error_temp))
)
write.csv(diagnostico_semanal,
          file.path(DIR_SALIDA, "diagnostico_temperatura_por_semana.csv"),
          row.names = FALSE)
write.csv(diagnostico_celdas,
          file.path(DIR_SALIDA, "diagnostico_temperatura_por_celda.csv"),
          row.names = FALSE)

# Cuantifica hipótesis, pero NO usa estas correcciones en los GeoTIFF candidatos.
ok <- is.finite(error_temp)
rmse_original <- sqrt(mean(error_temp[ok]^2))
error_sin_sesgo_celda <- sweep(error_temp, 1, sesgo_celda, "-")
error_sin_sesgo_semana <- sweep(error_temp, 2, sesgo_semana, "-")
resumen_diagnostico <- data.frame(
  escenario = c("bilineal_sin_correccion", "restar_sesgo_fijo_por_celda",
                "restar_sesgo_medio_por_semana"),
  RMSE = c(rmse_original,
           sqrt(mean(error_sin_sesgo_celda[is.finite(error_sin_sesgo_celda)]^2)),
           sqrt(mean(error_sin_sesgo_semana[is.finite(error_sin_sesgo_semana)]^2))),
  MAE = c(mean(abs(error_temp[ok])),
          mean(abs(error_sin_sesgo_celda[is.finite(error_sin_sesgo_celda)])),
          mean(abs(error_sin_sesgo_semana[is.finite(error_sin_sesgo_semana)])))
)
write.csv(resumen_diagnostico,
          file.path(DIR_SALIDA, "prueba_hipotesis_sesgo_temperatura.csv"),
          row.names = FALSE)

# Productos candidatos primarios y alternativa para sensibilidad.
writeRaster(candidatos_temp$bilinear,
            file.path(DIR_SALIDA, "temperatura_2024_bilineal_candidata.tif"),
            overwrite = TRUE, datatype = "FLT4S")
writeRaster(candidatos_rad$near,
            file.path(DIR_SALIDA, "radiacion_2024_vecino_candidata.tif"),
            overwrite = TRUE, datatype = "FLT4S")
writeRaster(candidatos_rad$bilinear,
            file.path(DIR_SALIDA, "radiacion_2024_bilineal_sensibilidad.tif"),
            overwrite = TRUE, datatype = "FLT4S")

mapa_sesgo <- plantilla
values(mapa_sesgo) <- NA_real_
mapa_sesgo[celdas_valle] <- sesgo_celda
writeRaster(mapa_sesgo,
            file.path(DIR_SALIDA, "sesgo_medio_temperatura_entregado_menos_nasa.tif"),
            overwrite = TRUE, datatype = "FLT4S")

cat("\n", strrep("=", 79), "\n", sep = "")
cat("RECONSTRUCCION REGIONAL NASA POWER 2024\n")
cat("52 bloques de 7 dias: 1 enero a 29 diciembre; LST; comunidad AG\n")
cat("Soporte aprobado:", length(celdas_valle), "celdas por centro de pixel\n")
cat(strrep("=", 79), "\n", sep = "")
print(tabla_metricas, row.names = FALSE, digits = 5)
cat("\nPRUEBA DEL PATRON DE ERROR DE TEMPERATURA\n")
print(resumen_diagnostico, row.names = FALSE, digits = 5)
cat("\nLectura: si quitar el sesgo fijo por celda reduce mucho el RMSE, la\n")
cat("diferencia es principalmente espacial y estable, no una semana mal armada.\n")
cat("Rango del sesgo semanal:", round(range(sesgo_semana), 4), "grados C\n")
cat("Rango del sesgo por celda:", round(range(sesgo_celda, na.rm = TRUE), 4),
    "grados C\n")

dibujar_resultados <- function() {
  # Pagina 1: comparación visual semana 10.
  par(mfrow = c(2, 3), mar = c(2.7, 2.7, 3, 4))
  plot(temp_entregada[[10]], main = "Temp. entregada S10", axes = FALSE)
  lines(valle)
  plot(candidatos_temp$bilinear[[10]], main = "Temp. bilineal S10", axes = FALSE)
  lines(valle)
  plot(temp_entregada[[10]] - candidatos_temp$bilinear[[10]],
       main = "Error temperatura S10", axes = FALSE)
  lines(valle)
  plot(rad_entregada[[10]], main = "Radiacion entregada S10", axes = FALSE)
  lines(valle)
  plot(candidatos_rad$near[[10]], main = "Radiacion vecino S10", axes = FALSE)
  lines(valle)
  plot(rad_entregada[[10]] - candidatos_rad$near[[10]],
       main = "Error radiacion S10", axes = FALSE)
  lines(valle)

  # Pagina 2: evidencia de la decisión entre métodos.
  par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))
  for (variable in c("temperatura", "radiacion")) {
    z <- tabla_metricas[tabla_metricas$variable == variable &
                         tabla_metricas$periodo == "semanas_02_52", ]
    barplot(z$RMSE, names.arg = z$metodo, col = "steelblue",
            ylab = "RMSE", main = paste("RMSE", variable, "S02-S52"))
  }
  plot(1:52, sesgo_semana, type = "b", pch = 19, col = "firebrick",
       xlab = "Semana", ylab = "Entregado - NASA (°C)",
       main = "Sesgo medio semanal temperatura")
  abline(h = 0, lty = 2)
  plot(mapa_sesgo, main = "Sesgo medio por celda (°C)", axes = FALSE)
  lines(valle)

  # Pagina 3: temporalidad y estabilidad espacial.
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(1:52, sd_espacial_semana, type = "b", pch = 19,
       xlab = "Semana", ylab = "SD espacial del error (°C)",
       main = "Variación espacial del error")
  hist(sesgo_celda, breaks = 20, col = "tan", border = "white",
       xlab = "Sesgo medio (°C)", main = "Sesgo fijo por celda")
  boxplot(error_temp[, 2:52], outline = FALSE, col = "lightblue",
          xaxt = "n", xlab = "Semanas 2-52", ylab = "Error (°C)",
          main = "Error de temperatura por semana")
  axis(1, at = seq(1, 51, 10), labels = seq(2, 52, 10))
  barplot(resumen_diagnostico$RMSE,
          names.arg = c("Original", "Sin sesgo\ncelda", "Sin sesgo\nsemana"),
          col = c("firebrick", "darkgreen", "goldenrod"), ylab = "RMSE (°C)",
          main = "Prueba diagnóstica (no corrección)")
}

# Se muestra en RStudio Plots y se conserva una copia multipágina.
dibujar_resultados()
pdf(file.path(DIR_SALIDA, "diagnostico_reconstruccion_power_2024.pdf"),
    width = 11, height = 8.5, onefile = TRUE)
dibujar_resultados()
dev.off()

cat("\nCandidatos completos y diagnosticos guardados en:\n", DIR_SALIDA, "\n")
cat("IMPORTANTE: no se aplico ninguna correccion de sesgo a temperatura.\n")
cat("Decision pendiente antes de 2010-2025: aprobar temperatura bilineal,\n")
cat("radiacion por vecino como principal y bilineal como sensibilidad.\n")
