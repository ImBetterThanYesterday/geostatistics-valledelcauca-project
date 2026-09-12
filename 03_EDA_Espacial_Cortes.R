# =============================================================================
# 03_EDA_Espacial_Cortes.R
#
# FASES 1 Y 2 DEL ANALISIS ESPACIAL
#
# Fase 1: control de calidad de los tres cortes seleccionados.
# Fase 2: EDA descriptivo espacial mediante mapas y relaciones simples.
#
# Todavia NO se ajustan modelos, no se calculan residuales y no se construyen
# semivariogramas. Esos pasos comienzan despues de revisar estos resultados.
# =============================================================================

suppressPackageStartupMessages(library(terra))

buscar_base <- function() {
  candidatos <- unique(normalizePath(c(
    getwd(), file.path(getwd(), ".."), file.path(getwd(), "../.."),
    file.path(getwd(), "../../..")
  ), mustWork = FALSE))
  ok <- vapply(candidatos, function(p) {
    file.exists(file.path(p, "avance_team_proyecto1", "reinicio_desde_cero",
                          "cortes_espaciales_resultados",
                          "tres_cortes_espaciales.rds"))
  }, logical(1))
  if (!any(ok)) stop("No se encontró tres_cortes_espaciales.rds. Ejecute antes ",
                     "02_01_Seleccion_Cortes_Espaciales.R.")
  candidatos[which(ok)[1]]
}

BASE <- buscar_base()
DIR_CORTES <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                        "cortes_espaciales_resultados")
DIR_SALIDA <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                        "eda_espacial_resultados")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

cortes <- readRDS(file.path(DIR_CORTES, "tres_cortes_espaciales.rds"))
categorias <- c("seca", "intermedia", "humeda")
cortes <- cortes[categorias]

columnas <- c("id_celda", "x", "y", "anio", "semana", "altitud_m",
              "precipitacion_mm", "temperatura_bilineal_c",
              "radiacion_bilineal_mj_m2_dia")

# -----------------------------------------------------------------------------
# FASE 1: control de calidad por corte
# -----------------------------------------------------------------------------
control <- do.call(rbind, lapply(names(cortes), function(cat) {
  d <- cortes[[cat]]
  faltantes <- setdiff(columnas, names(d))
  if (length(faltantes) > 0) stop(cat, ": faltan columnas: ",
                                  paste(faltantes, collapse = ", "))
  data.frame(
    categoria = cat,
    anio = unique(d$anio),
    banda = unique(d$semana),
    n_filas = nrow(d),
    n_celdas = length(unique(d$id_celda)),
    duplicados_id = sum(duplicated(d$id_celda)),
    duplicados_xy = sum(duplicated(d[, c("x", "y")])),
    nulos = sum(is.na(d[, columnas])),
    precip_negativa = sum(d$precipitacion_mm < 0, na.rm = TRUE),
    temp_fuera_rango = sum(d$temperatura_bilineal_c < -10 |
                             d$temperatura_bilineal_c > 45, na.rm = TRUE),
    rad_fuera_rango = sum(d$radiacion_bilineal_mj_m2_dia < 0 |
                            d$radiacion_bilineal_mj_m2_dia > 45, na.rm = TRUE),
    altitud_min = min(d$altitud_m, na.rm = TRUE),
    altitud_max = max(d$altitud_m, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rownames(control) <- NULL

if (any(control$n_filas != 686 | control$n_celdas != 686 |
        control$duplicados_id > 0 | control$duplicados_xy > 0 |
        control$nulos > 0 | control$precip_negativa > 0 |
        control$temp_fuera_rango > 0 | control$rad_fuera_rango > 0)) {
  stop("El control de calidad detectó problemas. Revise control_cortes.csv.")
}

write.csv(control, file.path(DIR_SALIDA, "control_cortes.csv"), row.names = FALSE)
cat("CONTROL DE CALIDAD\n")
print(control, row.names = FALSE)

# Resumen descriptivo de las variables en cada corte.
variables <- c("precipitacion_mm", "temperatura_bilineal_c",
               "radiacion_bilineal_mj_m2_dia", "altitud_m")
resumen <- do.call(rbind, lapply(names(cortes), function(cat) {
  d <- cortes[[cat]]
  do.call(rbind, lapply(variables, function(v) {
    x <- d[[v]]
    data.frame(categoria = cat, variable = v, n = sum(is.finite(x)),
               media = mean(x), mediana = median(x), sd = sd(x),
               minimo = min(x), maximo = max(x),
               p25 = unname(quantile(x, .25)),
               p75 = unname(quantile(x, .75)), stringsAsFactors = FALSE)
  }))
}))
rownames(resumen) <- NULL
write.csv(resumen, file.path(DIR_SALIDA, "resumen_variables_por_corte.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# FASE 2: EDA espacial y visualización
# -----------------------------------------------------------------------------
crear_raster <- function(d, variable) {
  terra::rast(d[, c("x", "y", variable)], type = "xyz", crs = "EPSG:4326")
}

pdf(file.path(DIR_SALIDA, "EDA_espacial_tres_cortes.pdf"),
    width = 11, height = 8.5)

for (cat in names(cortes)) {
  d <- cortes[[cat]]
  etiqueta <- paste0(toupper(substr(cat, 1, 1)), substr(cat, 2, nchar(cat)),
                     " | ", unique(d$anio), " - banda ", unique(d$semana))

  # Pagina 1: cuatro superficies espaciales.
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.2, 4.8))
  for (v in variables) {
    r <- crear_raster(d, v)
    titulo <- switch(v,
      precipitacion_mm = "Precipitación semanal (mm)",
      temperatura_bilineal_c = "Temperatura bilineal (°C)",
      radiacion_bilineal_mj_m2_dia = "Radiación bilineal (MJ/m²/día)",
      altitud_m = "Altitud (m)")
    pal <- switch(v,
      precipitacion_mm = hcl.colors(40, "Blues 3", rev = TRUE),
      temperatura_bilineal_c = hcl.colors(40, "Inferno"),
      radiacion_bilineal_mj_m2_dia = hcl.colors(40, "YlOrBr"),
      altitud_m = hcl.colors(40, "Terrain 2"))
    plot(r, main = titulo, xlab = "Longitud", ylab = "Latitud", col = pal)
  }
  mtext(paste("EDA espacial:", etiqueta), outer = TRUE, line = -1.2,
        font = 2, cex = 1.2)

  # Página 2: relaciones de precipitación con cada covariable.
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.2, 1.5))
  covariables <- c("altitud_m", "temperatura_bilineal_c",
                   "radiacion_bilineal_mj_m2_dia")
  colores <- c("darkgreen", "firebrick", "darkorange")
  for (j in seq_along(covariables)) {
    v <- covariables[j]
    plot(d[[v]], d$precipitacion_mm, pch = 19, cex = .55,
         col = adjustcolor(colores[j], .5),
         xlab = v, ylab = "Precipitación (mm)",
         main = paste("Precipitación vs", v))
    abline(lm(d$precipitacion_mm ~ d[[v]]), col = "black", lwd = 2)
    grid()
  }
  hist(d$precipitacion_mm, breaks = 30, col = "steelblue", border = "white",
       main = "Distribución de precipitación", xlab = "mm", ylab = "Frecuencia")
  mtext(paste("Relaciones y distribución:", etiqueta), outer = TRUE, line = -1.2,
        font = 2, cex = 1.2)
}
dev.off()

cat("\nFases 1 y 2 terminadas. Productos guardados en:\n", DIR_SALIDA, "\n", sep = "")
cat("- control_cortes.csv\n")
cat("- resumen_variables_por_corte.csv\n")
cat("- EDA_espacial_tres_cortes.pdf\n")
cat("Revise el PDF antes de pasar a modelos de media.\n")
