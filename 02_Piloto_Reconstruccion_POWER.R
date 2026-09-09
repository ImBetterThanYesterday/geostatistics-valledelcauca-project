# =============================================================================
# PROYECTO 1 - PILOTO FORENSE PARA RECONSTRUIR NASA POWER
#
# Objetivo: identificar, antes de descargar 2010-2025, la receta que produjo
# los valores semanales suministrados. No imputa ni modifica ningún GeoTIFF.
# Decisiones y resultados deben registrarse en BITACORA_PROYECTO.md.
# =============================================================================

library(terra)

ANIO_PILOTO <- 2024L
TOLERANCIA_COINCIDENCIA <- 1e-4

detectar_raiz_proyecto <- function() {
  inicios <- getwd()
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    activo <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                       error = function(e) "")
    if (nzchar(activo)) inicios <- c(inicios, dirname(activo))
  }
  archivo_source <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (length(archivo_source) && nzchar(archivo_source)) {
    inicios <- c(inicios, dirname(normalizePath(archivo_source)))
  }
  for (inicio in unique(inicios)) {
    actual <- normalizePath(inicio, mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(actual, "datos_proyecto_1"))) return(actual)
      padre <- dirname(actual)
      if (identical(padre, actual)) break
      actual <- padre
    }
  }
  stop("No se encontro el proyecto. Abra Proyecto1.Rproj.")
}

RAIZ <- detectar_raiz_proyecto()
setwd(RAIZ)
DIR_PILOTO <- file.path(RAIZ, "avance_team_proyecto1", "reinicio_desde_cero",
                        "power_piloto")
DIR_RAW <- file.path(DIR_PILOTO, "raw_nasa")
DIR_RESULTADOS <- file.path(DIR_PILOTO, "resultados")
dir.create(DIR_RAW, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_RESULTADOS, recursive = TRUE, showWarnings = FALSE)

# Nodos elegidos porque contienen valores suministrados de temperatura,
# radiacion o ambas. Los CSV diarios crudos se conservan sin editar.
puntos <- data.frame(
  id = sprintf("punto_%02d", 1:3),
  celda = c(1260L, 138L, 1094L),
  lon = c(-77.475, -76.225, -76.525),
  lat = c(3.275, 4.825, 3.525),
  archivo = sprintf("punto_%02d_2024.csv", 1:3)
)

crear_url <- function(lon, lat) {
  paste0(
    "https://power.larc.nasa.gov/api/temporal/daily/point?",
    "parameters=T2M,ALLSKY_SFC_SW_DWN&community=AG",
    "&longitude=", lon, "&latitude=", lat,
    "&start=20240101&end=20241231&format=CSV&time-standard=LST"
  )
}
puntos$url <- mapply(crear_url, puntos$lon, puntos$lat)
write.csv(puntos, file.path(DIR_RESULTADOS, "puntos_y_urls_piloto.csv"),
          row.names = FALSE)

# Descargar solo si el archivo crudo no existe. Esto evita solicitudes
# repetidas y mantiene la respuesta original para auditoria.
for (i in seq_len(nrow(puntos))) {
  destino <- file.path(DIR_RAW, puntos$archivo[i])
  if (!file.exists(destino)) {
    cat("Descargando", puntos$id[i], "desde NASA POWER...\n")
    tryCatch(download.file(puntos$url[i], destino, mode = "wb", quiet = FALSE),
             error = function(e) stop("No se pudo descargar ", puntos$id[i],
                                      ": ", conditionMessage(e)))
  } else {
    cat("Usando archivo crudo existente:", basename(destino), "\n")
  }
}

leer_power_csv <- function(archivo) {
  lineas <- readLines(archivo, warn = FALSE)
  fila <- grep("^YEAR,", lineas)
  if (length(fila) != 1) stop("Cabecera NASA no reconocida: ", archivo)
  datos <- read.csv(text = paste(lineas[fila:length(lineas)], collapse = "\n"))
  attr(datos, "cabecera") <- lineas[seq_len(fila - 1)]
  datos
}

agregar_52_semanas <- function(diario) {
  # Receta candidata: bloques consecutivos de 7 dias desde el 1 de enero.
  # En 2024 usa DOY 1-364 y deja fuera DOY 365-366.
  diario$semana_bloque <- (diario$DOY - 1L) %/% 7L + 1L
  usado <- diario[diario$semana_bloque <= 52L, ]
  salida <- aggregate(cbind(T2M, ALLSKY_SFC_SW_DWN) ~ semana_bloque,
                      data = usado, FUN = mean, na.rm = TRUE)
  names(salida) <- c("semana", "temperatura_nasa", "radiacion_nasa")
  salida
}

r_power <- rast(file.path(RAIZ, "datos_proyecto_1", "imagenes_semanales",
                           sprintf("power_semanal_valle_%d.tif", ANIO_PILOTO)))
if (nlyr(r_power) != 104L) stop("El archivo POWER no contiene 104 bandas.")

comparaciones <- data.frame()
metadatos <- data.frame()
for (i in seq_len(nrow(puntos))) {
  archivo <- file.path(DIR_RAW, puntos$archivo[i])
  diario <- leer_power_csv(archivo)
  semanal <- agregar_52_semanas(diario)
  temp_entregada <- as.numeric(values(r_power[[1:52]], mat = TRUE)[puntos$celda[i], ])
  rad_entregada <- as.numeric(values(r_power[[53:104]], mat = TRUE)[puntos$celda[i], ])
  comparaciones <- rbind(comparaciones, data.frame(
    id = puntos$id[i], celda = puntos$celda[i],
    lon = puntos$lon[i], lat = puntos$lat[i], semana = semanal$semana,
    temperatura_entregada = temp_entregada,
    temperatura_nasa = semanal$temperatura_nasa,
    radiacion_entregada = rad_entregada,
    radiacion_nasa = semanal$radiacion_nasa
  ))
  cab <- attr(diario, "cabecera")
  metadatos <- rbind(metadatos, data.frame(
    id = puntos$id[i],
    fuente = cab[grep("NASA/POWER Source", cab)][1],
    localizacion = cab[grep("^Location:", cab)][1],
    temperatura = cab[grep("^T2M ", cab)][1],
    radiacion = cab[grep("^ALLSKY_SFC_SW_DWN", cab)][1]
  ))
}
write.csv(comparaciones,
          file.path(DIR_RESULTADOS, "comparacion_semanal_2024.csv"),
          row.names = FALSE)
write.csv(metadatos, file.path(DIR_RESULTADOS, "metadatos_nasa.csv"),
          row.names = FALSE)

calcular_metricas <- function(entregado, nasa, variable, punto) {
  ok <- is.finite(entregado) & is.finite(nasa)
  if (!any(ok)) return(NULL)
  error <- entregado[ok] - nasa[ok]
  data.frame(
    punto = punto, variable = variable, n = sum(ok),
    ME = mean(error), MAE = mean(abs(error)), RMSE = sqrt(mean(error^2)),
    correlacion = if (sum(ok) > 2) cor(entregado[ok], nasa[ok]) else NA_real_,
    diferencia_maxima = max(abs(error)),
    coincidencias_exactas = sum(abs(error) <= TOLERANCIA_COINCIDENCIA),
    porcentaje_coincidencia = 100 * mean(abs(error) <= TOLERANCIA_COINCIDENCIA)
  )
}

metricas <- data.frame()
for (id in unique(comparaciones$id)) {
  z <- comparaciones[comparaciones$id == id, ]
  metricas <- rbind(metricas,
    calcular_metricas(z$temperatura_entregada, z$temperatura_nasa,
                      "temperatura", id),
    calcular_metricas(z$radiacion_entregada, z$radiacion_nasa,
                      "radiacion", id))
}
write.csv(metricas, file.path(DIR_RESULTADOS, "metricas_coincidencia.csv"),
          row.names = FALSE)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("PILOTO FORENSE NASA POWER —", ANIO_PILOTO, "\n")
cat("Receta candidata: AG + LST + bloques de 7 dias desde el 1 de enero\n")
cat("Temperatura: T2M, promedio semanal, grados C\n")
cat("Radiacion: ALLSKY_SFC_SW_DWN, promedio semanal, MJ/m2/dia\n")
cat(strrep("=", 78), "\n", sep = "")
print(metricas, row.names = FALSE)

semana_control <- comparaciones[comparaciones$semana == 10, ]
cat("\nCOMPARACION PARTICULAR DE LA SEMANA 10\n")
print(semana_control[, c("id", "temperatura_entregada", "temperatura_nasa",
                         "radiacion_entregada", "radiacion_nasa")],
      row.names = FALSE)

dibujar_piloto <- function() {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (variable in c("temperatura", "radiacion")) {
    col_e <- paste0(variable, "_entregada")
    col_n <- paste0(variable, "_nasa")
    z <- comparaciones[is.finite(comparaciones[[col_e]]), ]
    plot(z[[col_n]], z[[col_e]], pch = 19,
         col = if (variable == "temperatura") "firebrick" else "darkorange",
         xlab = paste(variable, "NASA reconstruida"),
         ylab = paste(variable, "GeoTIFF entregado"),
         main = paste("Coincidencia anual -", variable))
    abline(0, 1, lty = 2, col = "black")
  }
  for (variable in c("temperatura", "radiacion")) {
    col_e <- paste0(variable, "_entregada")
    col_n <- paste0(variable, "_nasa")
    z <- comparaciones[comparaciones$id == if (variable == "temperatura")
                         "punto_02" else "punto_03", ]
    matplot(z$semana, cbind(z[[col_e]], z[[col_n]]), type = "l", lwd = 2,
            lty = c(1, 2), col = c("navy", "red"),
            xlab = "Semana", ylab = variable,
            main = paste("Serie semanal -", variable))
    legend("topright", c("Entregado", "NASA actual"),
           col = c("navy", "red"), lty = c(1, 2), cex = .8)
  }

  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  for (variable in c("temperatura", "radiacion")) {
    col_e <- paste0(variable, "_entregada")
    col_n <- paste0(variable, "_nasa")
    z <- comparaciones[is.finite(comparaciones[[col_e]]), ]
    error <- z[[col_e]] - z[[col_n]]
    plot(z$semana, error, pch = 19, col = rgb(.2, .3, .7, .5),
         xlab = "Semana", ylab = "Entregado - NASA",
         main = paste("Error por semana -", variable))
    abline(h = 0, lty = 2, col = "red")
  }
}

# Visible en el panel Plots.
dibujar_piloto()

# Copia persistente.
pdf(file.path(DIR_RESULTADOS, "graficos_piloto_power.pdf"),
    width = 11, height = 8.5, onefile = TRUE)
dibujar_piloto()
dev.off()

cat("\nArchivos crudos preservados en:", DIR_RAW, "\n")
cat("Resultados guardados en:", DIR_RESULTADOS, "\n")
cat("El piloto no modifica ni completa los GeoTIFF suministrados.\n")
