# =============================================================================
# PROYECTO 1 - FASE 2: SOPORTE ESPACIAL-TEMPORAL COMUN
# Decisiones aprobadas:
#   D-001. Mascara oficial: centro del pixel dentro del Valle (686 celdas).
#   D-002. POWER no se imputa; se conserva su disponibilidad original.
#   D-003. Soporte objetivo independiente de disponibilidad de covariables.
# Trazabilidad completa: BITACORA_PROYECTO.md.
# =============================================================================

library(terra)

ANIO_REVISION <- 2024L
SEMANA_REVISION <- 10L

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
DIR_DATOS <- file.path(RAIZ, "datos_proyecto_1", "imagenes_semanales")
DIR_GADM <- file.path(RAIZ, "datos_proyecto_1", "cache_raster", "gadm",
                      "gadm41_COL_1_pk.rds")
DIR_SALIDA <- file.path(RAIZ, "avance_team_proyecto1", "reinicio_desde_cero",
                        "soporte_resultados")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ Archivos -------------------------------------
r_prec_anual <- rast(file.path(DIR_DATOS,
                                sprintf("chirps_semanal_valle_%d.tif", ANIO_REVISION)))
r_power_anual <- rast(file.path(DIR_DATOS,
                                 sprintf("power_semanal_valle_%d.tif", ANIO_REVISION)))
r_alt <- rast(file.path(DIR_DATOS, "altitud_valle.tif"))
r_clim_prec <- rast(file.path(DIR_DATOS, "chirps_climatologia_semanal_valle.tif"))
r_clim_temp <- rast(file.path(DIR_DATOS, "power_temp_climatologia_semanal_valle.tif"))
r_clim_rad <- rast(file.path(DIR_DATOS, "power_radiacion_climatologia_semanal_valle.tif"))

r_prec <- r_prec_anual[[SEMANA_REVISION]]
r_temp <- r_power_anual[[SEMANA_REVISION]]
r_rad <- r_power_anual[[52L + SEMANA_REVISION]]
names(r_prec) <- "precipitacion"
names(r_temp) <- "temperatura"
names(r_rad) <- "radiacion"
names(r_alt) <- "altitud"

gadm_obj <- readRDS(DIR_GADM)
colombia <- if (inherits(gadm_obj, "PackedSpatVector")) vect(gadm_obj) else gadm_obj
campo <- names(colombia)[grepl("NAME_1", names(colombia))][1]
valle <- colombia[grepl("Valle del Cauca", colombia[[campo]][[1]],
                        ignore.case = TRUE), ]

# Soporte maestro aprobado: centro de celda dentro del poligono GADM.
mascara_maestra <- rasterize(valle, r_prec, field = 1, background = NA,
                             touches = FALSE)
dentro <- !is.na(values(mascara_maestra, mat = FALSE))
celdas_maestras <- which(dentro)
xy <- xyFromCell(r_prec, celdas_maestras)

# Valores crudos: no se cambia ningun NA ni valor.
vp <- values(r_prec, mat = FALSE)[dentro]
vt <- values(r_temp, mat = FALSE)[dentro]
vr <- values(r_rad, mat = FALSE)[dentro]
va <- values(r_alt, mat = FALSE)[dentro]

# Los negativos CHIRPS son invalidos, pero dentro de la mascara centro no hay.
ok_p <- is.finite(vp) & vp >= 0
ok_t <- is.finite(vt)
ok_r <- is.finite(vr) & vr >= 0
ok_a <- is.finite(va)
ok_todos <- ok_p & ok_t & ok_r & ok_a

soporte <- data.frame(
  celda = celdas_maestras, x = xy[, 1], y = xy[, 2],
  dentro_valle = TRUE,
  precipitacion = vp, temperatura = vt, radiacion = vr, altitud = va,
  disponible_precipitacion = ok_p,
  disponible_temperatura = ok_t,
  disponible_radiacion = ok_r,
  disponible_altitud = ok_a,
  caso_completo_4_variables = ok_todos,
  anio = ANIO_REVISION, semana = SEMANA_REVISION
)
write.csv(soporte, file.path(DIR_SALIDA, "soporte_maestro_2024_S10.csv"),
          row.names = FALSE)
saveRDS(list(mascara = wrap(mascara_maestra), valle = wrap(valle),
             celdas = celdas_maestras, criterio = "centro del pixel"),
        file.path(DIR_SALIDA, "soporte_maestro.rds"))
writeRaster(mascara_maestra, file.path(DIR_SALIDA, "mascara_maestra_686.tif"),
            overwrite = TRUE, datatype = "INT1U")

# ------------------------ Diferencias entre datasets -------------------------
tabla_datasets <- data.frame(
  dataset = c("Precipitacion semanal", "Temperatura semanal",
              "Radiacion semanal", "Altitud",
              "Climatologia precipitacion", "Climatologia temperatura",
              "Climatologia radiacion"),
  fuente = c("CHIRPS", "NASA POWER", "NASA POWER", "SRTM",
             "CHIRPS", "NASA POWER", "NASA POWER"),
  archivo = c(basename(sources(r_prec_anual)[1]), basename(sources(r_power_anual)[1]),
              basename(sources(r_power_anual)[1]), basename(sources(r_alt)[1]),
              basename(sources(r_clim_prec)[1]), basename(sources(r_clim_temp)[1]),
              basename(sources(r_clim_rad)[1])),
  bandas_archivo = c(nlyr(r_prec_anual), nlyr(r_power_anual), nlyr(r_power_anual),
                     nlyr(r_alt), nlyr(r_clim_prec), nlyr(r_clim_temp), nlyr(r_clim_rad)),
  bandas_de_variable = c(52, 52, 52, 1, 52, 52, 52),
  soporte_temporal = c("semana-anio", "semana-anio", "semana-anio", "estatico",
                       "semana climatologica", "semana climatologica",
                       "semana climatologica"),
  periodo = c("2010-2025", "2010-2025", "2010-2025", "estatico",
              "promedio 2010-2025", "promedio 2010-2025", "promedio 2010-2025"),
  resolucion_grados = c(res(r_prec)[1], res(r_temp)[1], res(r_rad)[1], res(r_alt)[1],
                        res(r_clim_prec)[1], res(r_clim_temp)[1], res(r_clim_rad)[1]),
  pixeles_validos_semana = c(sum(ok_p), sum(ok_t), sum(ok_r), sum(ok_a),
                             sum(is.finite(values(r_clim_prec[[SEMANA_REVISION]], mat = FALSE)[dentro]) &
                                   values(r_clim_prec[[SEMANA_REVISION]], mat = FALSE)[dentro] >= 0),
                             sum(is.finite(values(r_clim_temp[[SEMANA_REVISION]], mat = FALSE)[dentro])),
                             sum(is.finite(values(r_clim_rad[[SEMANA_REVISION]], mat = FALSE)[dentro])))
)
tabla_datasets$cobertura_pct <- 100 * tabla_datasets$pixeles_validos_semana / length(celdas_maestras)
write.csv(tabla_datasets, file.path(DIR_SALIDA, "diferencias_datasets.csv"),
          row.names = FALSE)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("FASE 2 - SOPORTE ESPACIAL-TEMPORAL COMUN\n")
cat("Mascara aprobada: centro del pixel dentro del Valle\n")
cat("Soporte maestro:", length(celdas_maestras), "pixeles\n")
cat("Semana visualizada:", sprintf("%d-S%02d", ANIO_REVISION, SEMANA_REVISION), "\n")
cat(strrep("=", 78), "\n\n", sep = "")
print(tabla_datasets, row.names = FALSE)

cat("\nDISPONIBILIDAD SOBRE LOS 686 PIXELES\n")
tabla_disponibilidad <- data.frame(
  variable = c("precipitacion", "temperatura", "radiacion", "altitud",
               "interseccion de las cuatro"),
  validos = c(sum(ok_p), sum(ok_t), sum(ok_r), sum(ok_a), sum(ok_todos)),
  faltantes = c(sum(!ok_p), sum(!ok_t), sum(!ok_r), sum(!ok_a),
                length(ok_todos) - sum(ok_todos))
)
tabla_disponibilidad$cobertura_pct <-
  100 * tabla_disponibilidad$validos / length(celdas_maestras)
print(tabla_disponibilidad, row.names = FALSE)
write.csv(tabla_disponibilidad,
          file.path(DIR_SALIDA, "disponibilidad_soporte.csv"), row.names = FALSE)

# ------------------------------- Visuales ------------------------------------
mapa_disponibilidad <- function(ok, nombre) {
  z <- mascara_maestra
  values(z) <- NA_real_
  z[celdas_maestras] <- as.integer(ok)
  plot(z, col = c("#ef3b2c", "#31a354"), breaks = c(-.5, .5, 1.5),
       legend = FALSE, axes = FALSE, main = nombre)
  plot(valle, add = TRUE, border = "black", lwd = 1.2)
  legend("bottomleft", c("NA original", "Disponible"),
         fill = c("#ef3b2c", "#31a354"), cex = .7, bty = "n")
}

dibujar_soporte <- function() {
  # Pagina 1: valores crudos en exactamente la misma extension.
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))
  plot(mask(r_prec, mascara_maestra), main = "CHIRPS: precipitacion semanal")
  plot(valle, add = TRUE)
  plot(mask(r_temp, mascara_maestra), main = "POWER: temperatura semanal")
  plot(valle, add = TRUE)
  plot(mask(r_rad, mascara_maestra), main = "POWER: radiacion semanal")
  plot(valle, add = TRUE)
  plot(mask(r_alt, mascara_maestra), main = "SRTM: altitud estatica")
  plot(valle, add = TRUE)

  # Pagina 2: cobertura comparable, verde disponible y rojo NA interno.
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 1))
  mapa_disponibilidad(ok_p, paste("Precipitacion:", sum(ok_p), "/ 686"))
  mapa_disponibilidad(ok_t, paste("Temperatura:", sum(ok_t), "/ 686"))
  mapa_disponibilidad(ok_r, paste("Radiacion:", sum(ok_r), "/ 686"))
  mapa_disponibilidad(ok_a, paste("Altitud:", sum(ok_a), "/ 686"))

  # Pagina 3: resumen y consecuencia para casos completos.
  par(mfrow = c(1, 2), mar = c(7, 4, 4, 1))
  barplot(tabla_disponibilidad$validos[1:4],
          names.arg = tabla_disponibilidad$variable[1:4], las = 2,
          col = c("#3182bd", "#e6550d", "#756bb1", "#31a354"),
          ylim = c(0, 720), ylab = "Pixeles disponibles",
          main = "Cobertura en soporte maestro")
  abline(h = 686, lty = 2)
  mapa_disponibilidad(ok_todos,
                      paste("Casos completos sin imputar:", sum(ok_todos), "/ 686"))

  # Pagina 4: diagrama conceptual de diferencias.
  par(mfrow = c(1, 1), mar = c(1, 1, 3, 1))
  plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 10))
  rect(0.4, 6.1, 3.0, 9.0, col = "#deebf7", border = "#3182bd", lwd = 2)
  rect(3.7, 6.1, 6.3, 9.0, col = "#fee6ce", border = "#e6550d", lwd = 2)
  rect(7.0, 6.1, 9.6, 9.0, col = "#e5f5e0", border = "#31a354", lwd = 2)
  text(1.7, 8.5, "CHIRPS", font = 2, cex = 1.2)
  text(1.7, 7.5, "Precipitacion\n52 semanas x 16 anos\n686/686 pixeles")
  text(5.0, 8.5, "NASA POWER", font = 2, cex = 1.2)
  text(5.0, 7.35, "Temperatura: 122/686\nRadiacion: 27/686\n52 semanas x 16 anos")
  text(8.3, 8.5, "SRTM", font = 2, cex = 1.2)
  text(8.3, 7.5, "Altitud\n1 capa estatica\n686/686 pixeles")
  arrows(1.7, 5.7, 4.3, 3.9, length = .1, lwd = 2)
  arrows(5.0, 5.7, 5.0, 3.9, length = .1, lwd = 2)
  arrows(8.3, 5.7, 5.7, 3.9, length = .1, lwd = 2)
  rect(3.1, 1.0, 6.9, 3.8, col = "#f0f0f0", border = "black", lwd = 2)
  text(5, 3.25, "SOPORTE OBJETIVO", font = 2, cex = 1.2)
  text(5, 2.25, "686 pixeles x semana", cex = 1.1)
  text(5, 1.45, paste("Interseccion observada sin imputar:", sum(ok_todos), "pixeles"),
       col = "#cb181d", font = 2)
  title("Diferencias y relacion entre los datasets")
}

# Visible en el panel Plots.
dibujar_soporte()

# Copia persistente para revision.
pdf(file.path(DIR_SALIDA, "visual_diferencias_y_soporte.pdf"),
    width = 11, height = 8.5, onefile = TRUE)
dibujar_soporte()
dev.off()

cat("\nDECISION REGISTRADA\n")
cat("  El soporte objetivo queda fijado en 686 pixeles por semana.\n")
cat("  POWER permanece como NA donde no fue observado.\n")
cat("  No se construye aun una tabla de modelacion con cuatro variables,\n")
cat("  porque la interseccion cruda contiene solo", sum(ok_todos), "pixeles.\n")
cat("\nVisual guardada en:\n",
    file.path(DIR_SALIDA, "visual_diferencias_y_soporte.pdf"), "\n")
