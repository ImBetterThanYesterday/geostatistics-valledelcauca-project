# =============================================================================
# PROYECTO 1 - REINICIO DESDE CERO
# FASE 1: auditoria de nulos y control de calidad (sin imputar ni borrar)
# Trazabilidad: BITACORA_PROYECTO.md, hechos H-001 a H-008 y D-001/D-002.
#
# Salidas simultaneas:
#   - consola: tablas y hallazgos;
#   - panel Plots de RStudio: todos los graficos quedan en el historial;
#   - carpeta auditoria_resultados: CSV, PDF y GeoTIFF de categorias.
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
  stop("No se encontro datos_proyecto_1. Abra Proyecto1.Rproj y vuelva a ejecutar.")
}

RAIZ <- detectar_raiz_proyecto()
setwd(RAIZ)
DIR_DATOS <- file.path(RAIZ, "datos_proyecto_1", "imagenes_semanales")
DIR_GADM <- file.path(RAIZ, "datos_proyecto_1", "cache_raster", "gadm",
                      "gadm41_COL_1_pk.rds")
DIR_SALIDA <- file.path(RAIZ, "avance_team_proyecto1", "reinicio_desde_cero",
                        "auditoria_resultados")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("FASE 1 - AUDITORIA SIN MODIFICAR LOS DATOS\n")
cat("Raiz:", RAIZ, "\n")
cat("Caso visual:", sprintf("%d-S%02d", ANIO_REVISION, SEMANA_REVISION), "\n")
cat("Salidas:", DIR_SALIDA, "\n")
cat(strrep("=", 78), "\n\n", sep = "")

# -------------------------- Lectura de datos crudos ---------------------------
archivo_prec <- file.path(DIR_DATOS,
                           sprintf("chirps_semanal_valle_%d.tif", ANIO_REVISION))
archivo_power <- file.path(DIR_DATOS,
                            sprintf("power_semanal_valle_%d.tif", ANIO_REVISION))
r_prec <- rast(archivo_prec)[[SEMANA_REVISION]]
r_power <- rast(archivo_power)
r_temp <- r_power[[SEMANA_REVISION]]
r_rad <- r_power[[52L + SEMANA_REVISION]]
r_alt <- rast(file.path(DIR_DATOS, "altitud_valle.tif"))
names(r_prec) <- "precipitacion"
names(r_temp) <- "temperatura"
names(r_rad) <- "radiacion"
names(r_alt) <- "altitud"

gadm_obj <- readRDS(DIR_GADM)
colombia <- if (inherits(gadm_obj, "PackedSpatVector")) vect(gadm_obj) else gadm_obj
campo <- names(colombia)[grepl("NAME_1", names(colombia))][1]
valle <- colombia[grepl("Valle del Cauca", colombia[[campo]][[1]],
                        ignore.case = TRUE), ]

# Dos definiciones candidatas del soporte. No se escoge ninguna todavia.
mascara_centro <- rasterize(valle, r_prec, field = 1, background = NA,
                            touches = FALSE)
mascara_interseccion <- rasterize(valle, r_prec, field = 1, background = NA,
                                  touches = TRUE)
# Proporcion aproximada del area de cada pixel que cae dentro del poligono.
fraccion_dentro <- rasterize(valle, r_prec, field = 1, background = 0,
                             cover = TRUE)
inside_centro <- !is.na(values(mascara_centro, mat = FALSE))
inside_interseccion <- !is.na(values(mascara_interseccion, mat = FALSE))
v_fraccion <- values(fraccion_dentro, mat = FALSE)
frontera_sin_centro <- inside_interseccion & !inside_centro

cat("COMPARACION DE MASCARAS\n")
cat("  Centro del pixel dentro del Valle :", sum(inside_centro), "pixeles\n")
cat("  Pixel intersecta/toca el Valle    :", sum(inside_interseccion), "pixeles\n")
cat("  Frontera adicional                :",
    sum(inside_interseccion & !inside_centro), "pixeles\n\n")
cat("De los pixeles fronterizos sin centro interior:\n")
cat("  Mediana del area dentro del Valle :",
    sprintf("%.1f%%\n", 100 * median(v_fraccion[frontera_sin_centro])))
cat("  Menos de 10% dentro               :",
    sum(v_fraccion[frontera_sin_centro] < .10), "pixeles\n")
cat("  Al menos 50% dentro               :",
    sum(v_fraccion[frontera_sin_centro] >= .50), "pixeles\n\n")

xy_interseccion <- xyFromCell(r_prec, which(inside_interseccion))
write.csv(data.frame(
  celda = which(inside_interseccion), x = xy_interseccion[, 1],
  y = xy_interseccion[, 2], fraccion_dentro = v_fraccion[inside_interseccion],
  centro_dentro = inside_centro[inside_interseccion]
), file.path(DIR_SALIDA, "fraccion_area_dentro_valle.csv"), row.names = FALSE)

# ------------------------ Clasificacion conservadora --------------------------
# Codigos estables durante todo el proyecto.
etiquetas <- c("Valido", "NA_estructural", "NA_interno_original",
               "NA_control_calidad", "Dato_fuera_Valle")
colores <- c("#2ca25f", "#eeeeee", "#fdae6b", "#de2d26", "#756bb1")

es_invalido_conocido <- function(x, variable) {
  # Solo reglas inequívocas en esta etapa. No se declara invalida la altitud
  # -2.9 m ni se imponen umbrales nuevos sin decision del equipo.
  switch(variable,
         precipitacion = is.finite(x) & x < 0,
         radiacion = is.finite(x) & x < 0,
         temperatura = rep(FALSE, length(x)),
         altitud = rep(FALSE, length(x)))
}

clasificar <- function(r, inside, variable) {
  x <- values(r, mat = FALSE)
  invalido <- es_invalido_conocido(x, variable)
  codigo <- integer(length(x))
  codigo[!inside & !is.finite(x)] <- 2L                  # NA estructural
  codigo[inside & !is.finite(x)] <- 3L                   # NA interno original
  codigo[inside & is.finite(x) & invalido] <- 4L         # control de calidad
  codigo[!inside & is.finite(x)] <- 5L                   # dato fuera del Valle
  codigo[inside & is.finite(x) & !invalido] <- 1L        # valido
  salida <- r
  values(salida) <- codigo
  levels(salida) <- data.frame(ID = 1:5, categoria = etiquetas)
  names(salida) <- paste0("categoria_", variable)
  salida
}

rasters <- list(precipitacion = r_prec, temperatura = r_temp,
                radiacion = r_rad, altitud = r_alt)
mascaras <- list(centro = inside_centro, interseccion = inside_interseccion)
clasificados <- list()
conteos <- data.frame()

for (criterio in names(mascaras)) {
  clasificados[[criterio]] <- list()
  for (variable in names(rasters)) {
    rc <- clasificar(rasters[[variable]], mascaras[[criterio]], variable)
    clasificados[[criterio]][[variable]] <- rc
    tab <- tabulate(values(rc, mat = FALSE), nbins = 5)
    conteos <- rbind(conteos,
                     data.frame(criterio = criterio, variable = variable,
                                categoria = etiquetas, n = tab))
    writeRaster(rc,
                file.path(DIR_SALIDA,
                          paste0("categorias_", criterio, "_", variable, ".tif")),
                overwrite = TRUE, datatype = "INT1U")
  }
}

write.csv(conteos, file.path(DIR_SALIDA, "conteos_categorias_semana.csv"),
          row.names = FALSE)

cat("CLASIFICACION CONSERVADORA - SEMANA DE REVISION\n")
for (criterio in names(mascaras)) {
  cat("\n", toupper(criterio), "\n", strrep("-", 72), "\n", sep = "")
  tabla <- xtabs(n ~ variable + categoria,
                 conteos[conteos$criterio == criterio, ])
  print(tabla[, etiquetas, drop = FALSE])
}

# --------------------- Auditoria temporal de cobertura ------------------------
cat("\nAUDITORIA TEMPORAL 2010-2025 (puede tardar unos segundos)\n")
cobertura <- data.frame()
for (anio in 2010:2025) {
  rp <- rast(file.path(DIR_DATOS, sprintf("chirps_semanal_valle_%d.tif", anio)))
  rw <- rast(file.path(DIR_DATOS, sprintf("power_semanal_valle_%d.tif", anio)))
  grupos <- list(precipitacion = rp, temperatura = rw[[1:52]],
                 radiacion = rw[[53:104]])
  for (variable in names(grupos)) {
    matriz <- values(grupos[[variable]], mat = TRUE)
    invalidos <- apply(matriz, 2, es_invalido_conocido, variable = variable)
    if (is.null(dim(invalidos))) invalidos <- matrix(invalidos, ncol = 1)
    for (criterio in names(mascaras)) {
      dentro <- mascaras[[criterio]]
      cobertura <- rbind(cobertura, data.frame(
        anio = anio, semana = 1:52, variable = variable, criterio = criterio,
        validos = colSums(is.finite(matriz) & !invalidos & dentro),
        NA_interno_original = colSums(!is.finite(matriz) & dentro),
        NA_control_calidad = colSums(is.finite(matriz) & invalidos & dentro),
        Dato_fuera_Valle = colSums(is.finite(matriz) & !dentro)
      ))
    }
  }
  cat("  Revisado", anio, "\n")
}
write.csv(cobertura, file.path(DIR_SALIDA, "cobertura_2010_2025.csv"),
          row.names = FALSE)

resumen_cobertura <- aggregate(
  cbind(validos, NA_interno_original, NA_control_calidad, Dato_fuera_Valle) ~
    variable + criterio,
  data = cobertura,
  FUN = function(x) c(min = min(x), mediana = median(x), max = max(x))
)
write.csv(resumen_cobertura,
          file.path(DIR_SALIDA, "resumen_cobertura_2010_2025.csv"),
          row.names = FALSE)
cat("\nRESUMEN TEMPORAL: MINIMO / MEDIANA / MAXIMO POR SEMANA\n")
print(resumen_cobertura)

# ----------------------------- Graficos --------------------------------------
dibujar_auditoria <- function() {
  par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))
  plot(mascara_centro, col = "#9ecae1", legend = FALSE,
       main = paste("Mascara por centro\n", sum(inside_centro), "pixeles"))
  plot(valle, add = TRUE, border = "black", lwd = 1.5)
  plot(mascara_interseccion, col = "#74c476", legend = FALSE,
       main = paste("Mascara por interseccion\n", sum(inside_interseccion), "pixeles"))
  plot(valle, add = TRUE, border = "black", lwd = 1.5)
  plot(fraccion_dentro * 100, col = hcl.colors(20, "YlGnBu"),
       main = "Porcentaje del pixel dentro del Valle")
  plot(valle, add = TRUE, border = "black", lwd = 1.5)
  r_frontera <- r_prec
  values(r_frontera) <- NA_real_
  r_frontera[which(frontera_sin_centro)] <-
    100 * v_fraccion[frontera_sin_centro]
  plot(r_frontera, col = hcl.colors(20, "Inferno"),
       main = "103 pixeles tocados con centro afuera")
  plot(valle, add = TRUE, border = "black", lwd = 1.5)

  for (criterio in names(clasificados)) {
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 5))
    for (variable in names(rasters)) {
      mapa_numerico <- clasificados[[criterio]][[variable]]
      levels(mapa_numerico) <- NULL
      plot(mapa_numerico, col = colores, breaks = seq(0.5, 5.5, by = 1),
           legend = FALSE, axes = FALSE,
           main = paste(variable, "-", criterio))
      plot(valle, add = TRUE, border = "black", lwd = 1)
      presentes <- sort(unique(values(mapa_numerico, mat = FALSE)))
      legend("topright", legend = etiquetas[presentes], fill = colores[presentes],
             cex = .55, bty = "n")
    }
  }

  for (criterio in names(mascaras)) {
    sub <- conteos[conteos$criterio == criterio &
                     conteos$categoria %in% c("Valido", "NA_interno_original",
                                               "NA_control_calidad",
                                               "Dato_fuera_Valle"), ]
    mat <- xtabs(n ~ categoria + variable, sub)
    par(mfrow = c(1, 1), mar = c(6, 4, 4, 2))
    barplot(mat, beside = FALSE,
            col = colores[match(rownames(mat), etiquetas)], las = 2,
            ylab = "Numero de pixeles",
            main = paste("Categorias de calidad - criterio", criterio))
    legend("topright", legend = rownames(mat),
           fill = colores[match(rownames(mat), etiquetas)], cex = .8)
  }

  centro <- cobertura[cobertura$criterio == "centro", ]
  par(mfrow = c(3, 1), mar = c(3, 4, 3, 1))
  for (variable in c("precipitacion", "temperatura", "radiacion")) {
    z <- centro[centro$variable == variable, ]
    tiempo <- z$anio + (z$semana - 1) / 52
    plot(tiempo, z$validos, type = "l", col = "#238b45", lwd = 1,
         xlab = "Anio", ylab = "Pixeles validos",
         main = paste("Cobertura interior por semana -", variable))
  }
  par(mfrow = c(1, 1))
}

# Primera llamada: se ve en el panel Plots y queda en su historial.
dibujar_auditoria()

# Segunda llamada: conserva una copia multipagina para revision posterior.
pdf(file.path(DIR_SALIDA, "graficos_auditoria_nulos.pdf"),
    width = 11, height = 8.5, onefile = TRUE)
dibujar_auditoria()
dev.off()

cat("\n", strrep("=", 78), "\n", sep = "")
cat("AUDITORIA TERMINADA SIN IMPUTAR NI ELIMINAR DATOS\n")
cat("Revise los graficos con las flechas del panel Plots de RStudio.\n")
cat("PDF:", file.path(DIR_SALIDA, "graficos_auditoria_nulos.pdf"), "\n")
cat("Conteos:", file.path(DIR_SALIDA, "conteos_categorias_semana.csv"), "\n")
cat("Cobertura:", file.path(DIR_SALIDA, "cobertura_2010_2025.csv"), "\n")
cat(strrep("=", 78), "\n", sep = "")
