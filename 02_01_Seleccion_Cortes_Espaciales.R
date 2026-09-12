# =============================================================================
# 02_01_Seleccion_Cortes_Espaciales.R
#
# Selección reproducible de tres superficies espaciales semanales.
#
# Objetivo:
#   El dataset contiene 686 celdas x 16 años x 52 bandas semanales. Como el
#   siguiente análisis será solamente espacial, no se promedian las 832 bandas
#   entre sí. Se resume cada banda únicamente para clasificarla como seca,
#   intermedia o húmeda, y después se conservan las 686 celdas del corte real.
#
# Regla predefinida:
#   - percentil 10 de la precipitación media del Valle: corte seco;
#   - percentil 50: corte intermedio;
#   - percentil 90: corte húmedo.
#
# El promedio regional se usa SOLO para seleccionar cortes. La modelación
# espacial posterior utiliza la precipitación de cada celda, sin promediarla
# entre años o semanas.
# =============================================================================

buscar_base <- function() {
  candidatos <- unique(normalizePath(c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "../.."),
    file.path(getwd(), "../../..")
  ), mustWork = FALSE))

  ok <- vapply(candidatos, function(p) {
    file.exists(file.path(p, "dataset_eda", "dataset_final_7_variables.rds")) ||
      file.exists(file.path(p, "avance_team_proyecto1", "reinicio_desde_cero",
                            "dataset_eda", "dataset_final_7_variables.rds"))
  }, logical(1))

  if (!any(ok)) {
    stop("No se encontró dataset_eda/dataset_final_7_variables.rds. ",
         "Abra el proyecto desde la carpeta Proyecto1 o ajuste la ruta.")
  }
  candidatos[which(ok)[1]]
}

BASE <- buscar_base()
RUTA_DATASET <- file.path(BASE, "dataset_eda", "dataset_final_7_variables.rds")
if (!file.exists(RUTA_DATASET)) {
  RUTA_DATASET <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                            "dataset_eda", "dataset_final_7_variables.rds")
}
DIR_SALIDA <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                        "cortes_espaciales_resultados")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

dataset <- readRDS(RUTA_DATASET)

columnas_requeridas <- c("id_celda", "x", "y", "anio", "semana",
                         "precipitacion_mm")
faltantes <- setdiff(columnas_requeridas, names(dataset))
if (length(faltantes) > 0) {
  stop("Faltan columnas requeridas: ", paste(faltantes, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 1. Control de estructura
# -----------------------------------------------------------------------------
conteo_celdas <- length(unique(dataset$id_celda))
conteo_anios <- length(unique(dataset$anio))
conteo_bandas <- length(unique(dataset$semana))

cat("Dataset cargado:", nrow(dataset), "filas x", ncol(dataset), "columnas\n")
cat("Celdas:", conteo_celdas, "| Años:", conteo_anios,
    "| Bandas por año:", conteo_bandas, "\n")

if (conteo_celdas != 686 || conteo_anios != 16 || conteo_bandas != 52) {
  warning("La estructura no coincide exactamente con 686 x 16 x 52. " ,
          "Revisar antes de continuar.")
}

# Cada año-banda debe tener una observación por celda.
control_cortes <- aggregate(
  id_celda ~ anio + semana, data = dataset,
  FUN = function(z) length(unique(z))
)
names(control_cortes)[3] <- "n_celdas"
if (any(control_cortes$n_celdas != conteo_celdas)) {
  stop("Hay cortes año-banda con un número de celdas distinto de ",
       conteo_celdas, ". No se deben seleccionar hasta corregirlos.")
}

# -----------------------------------------------------------------------------
# 2. Resumen regional de cada superficie (solo para clasificar)
# -----------------------------------------------------------------------------
resumen_832 <- aggregate(precipitacion_mm ~ anio + semana, data = dataset,
                         FUN = function(z) mean(z, na.rm = TRUE))
names(resumen_832)[3] <- "media"
agregar_estadistica <- function(nombre, fun) {
  out <- aggregate(precipitacion_mm ~ anio + semana, data = dataset, FUN = fun)
  names(out)[3] <- nombre
  out
}
resumen_832 <- Reduce(function(x, y) merge(x, y, by = c("anio", "semana")),
                      list(resumen_832,
                           agregar_estadistica("mediana", median),
                           agregar_estadistica("sd", sd),
                           agregar_estadistica("minimo", min),
                           agregar_estadistica("maximo", max),
                           agregar_estadistica("n_validos", function(z) sum(is.finite(z)))))
resumen_832$corte <- paste0(resumen_832$anio, "_banda_",
                            sprintf("%02d", resumen_832$semana))
resumen_832 <- resumen_832[order(resumen_832$media), ]
rownames(resumen_832) <- NULL

if (any(resumen_832$n_validos != conteo_celdas)) {
  stop("Al menos un corte tiene valores faltantes de precipitación.")
}

# -----------------------------------------------------------------------------
# 3. Selección de los cortes más cercanos a P10, P50 y P90
# -----------------------------------------------------------------------------
probabilidades <- c(seca = 0.10, intermedia = 0.50, humeda = 0.90)
objetivos <- as.numeric(quantile(resumen_832$media,
                                 probs = probabilidades, type = 7,
                                 names = FALSE))

seleccion <- do.call(rbind, lapply(seq_along(probabilidades), function(i) {
  distancia <- abs(resumen_832$media - objetivos[i])
  fila <- which(distancia == min(distancia))[1]
  data.frame(
    categoria = names(probabilidades)[i],
    percentil = unname(probabilidades[i]),
    objetivo_mm = objetivos[i],
    anio = resumen_832$anio[fila],
    semana = resumen_832$semana[fila],
    corte = resumen_832$corte[fila],
    media_regional_mm = resumen_832$media[fila],
    diferencia_objetivo_mm = resumen_832$media[fila] - objetivos[i],
    sd_espacial_mm = resumen_832$sd[fila],
    minimo_mm = resumen_832$minimo[fila],
    maximo_mm = resumen_832$maximo[fila],
    stringsAsFactors = FALSE
  )
}))
rownames(seleccion) <- NULL

cat("\nCortes seleccionados:\n")
print(seleccion, row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Guardar productos reproducibles
# -----------------------------------------------------------------------------
write.csv(resumen_832,
          file.path(DIR_SALIDA, "resumen_832_cortes.csv"),
          row.names = FALSE)
write.csv(seleccion,
          file.path(DIR_SALIDA, "cortes_seleccionados_P10_P50_P90.csv"),
          row.names = FALSE)

# Guardar los tres data frames completos: cada uno es un dataset espacial de
# 686 celdas y un único año-banda.
cortes <- lapply(seq_len(nrow(seleccion)), function(i) {
  s <- seleccion[i, ]
  sub <- dataset[dataset$anio == s$anio & dataset$semana == s$semana, ]
  sub <- sub[order(sub$id_celda), ]
  attr(sub, "categoria_corte") <- s$categoria
  attr(sub, "regla_seleccion") <- paste0("Percentil ", 100 * s$percentil,
                                           " de la media regional")
  sub
})
names(cortes) <- seleccion$categoria
saveRDS(cortes, file.path(DIR_SALIDA, "tres_cortes_espaciales.rds"))

# -----------------------------------------------------------------------------
# 5. Visualización: distribución y ubicación de los cortes seleccionados
# -----------------------------------------------------------------------------
pdf(file.path(DIR_SALIDA, "seleccion_cortes_espaciales.pdf"),
    width = 11, height = 8.5)

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))
hist(resumen_832$media, breaks = 35, col = "grey80", border = "white",
     main = "Precipitación media por corte",
     xlab = "Media de precipitación en las 686 celdas (mm)",
     ylab = "Número de cortes")
colores <- c(seca = "steelblue", intermedia = "darkgreen", humeda = "firebrick")
for (i in seq_len(nrow(seleccion))) {
  abline(v = seleccion$media_regional_mm[i], col = colores[seleccion$categoria[i]],
         lwd = 3, lty = 2)
}
legend("topright", legend = paste0(seleccion$categoria, ": ",
                                   seleccion$anio, "-b", seleccion$semana),
       col = colores[seleccion$categoria], lwd = 3, lty = 2, bty = "n")

plot(seq_len(nrow(resumen_832)), resumen_832$media, pch = 16, cex = 0.55,
     col = "grey55", xaxt = "n",
     xlab = "Cortes ordenados de menor a mayor precipitación",
     ylab = "Precipitación media regional (mm)",
     main = "Regla P10–P50–P90")
for (i in seq_len(nrow(seleccion))) {
  fila <- which(resumen_832$corte == seleccion$corte[i])
  points(fila, resumen_832$media[fila], pch = 19, cex = 1.4,
         col = colores[seleccion$categoria[i]])
  text(fila, resumen_832$media[fila], labels = paste0(seleccion$categoria[i],
       "\n", seleccion$anio[i], "-b", seleccion$semana[i]),
       pos = 3, cex = 0.75, col = colores[seleccion$categoria[i]])
}
axis(1, at = c(1, nrow(resumen_832)), labels = c("más seco", "más húmedo"))

dev.off()

cat("\nProductos guardados en:\n", DIR_SALIDA, "\n", sep = "")
cat("- resumen_832_cortes.csv\n")
cat("- cortes_seleccionados_P10_P50_P90.csv\n")
cat("- tres_cortes_espaciales.rds\n")
cat("- seleccion_cortes_espaciales.pdf\n")
