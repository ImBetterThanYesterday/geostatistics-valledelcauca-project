###################################################################
# EDA por variable - dataset final
# Visualización en RStudio, sin guardar archivos
###################################################################

library(terra)

dataset <- readRDS("dataset_eda/dataset_final_7_variables.rds")
str(dataset)

variables <- c(
  "precipitacion_mm", "temperatura_bilineal_c",
  "radiacion_bilineal_mj_m2_dia", "altitud_m",
  "clim_precipitacion_mm", "clim_temperatura_bilineal_c",
  "clim_radiacion_bilineal_mj_m2_dia"
)

nombres <- c(
  "Precipitación", "Temperatura", "Radiación solar", "Altitud",
  "Climatología de precipitación", "Climatología de temperatura",
  "Climatología de radiación"
)

unidades <- c(
  "mm por periodo", "°C", "MJ/m²/día", "m",
  "mm por periodo", "°C", "MJ/m²/día"
)

nombres_matriz <- c(
  "Precipitación", "Temperatura", "Radiación", "Altitud",
  "Clim. precip.", "Clim. temp.", "Clim. rad."
)

colores <- c(
  "steelblue", "firebrick", "orange", "darkgreen",
  "steelblue4", "firebrick4", "darkorange3"
)

################################################################
# 0) Cobertura temporal aproximada
################################################################

# Esta aproximación no verifica las fechas reales de las bandas.
dataset$mes_aprox <- pmin(12, ceiling(dataset$semana * 7 / 30.4))
tabla_meses <- table(dataset$anio, dataset$mes_aprox) > 0

print(tabla_meses)
print(rowSums(tabla_meses))

################################################################
# 1) Estadísticas descriptivas y revisión de calidad
################################################################

# Rangos orientativos: estar fuera señala un valor para revisar.
rango_min <- c(0, -10, 0, -50, 0, -10, 0)
rango_max <- c(700, 45, 45, 6000, 700, 45, 45)

for (i in seq_along(variables)) {
  v <- variables[i]
  x <- dataset[[v]]
  
  cat("\n====================================================\n")
  cat(nombres[i], "|", unidades[i], "\n")
  cat("====================================================\n")
  print(summary(x))
  
  cat("Desviación estándar:", sd(x, na.rm = TRUE), "\n")
  cat("Nulos:", sum(is.na(x)), "\n")
  cat("Negativos:", sum(x < 0, na.rm = TRUE), "\n")
  cat("Ceros:", sum(x == 0, na.rm = TRUE),
      sprintf("(%.2f%%)\n", 100 * mean(x == 0, na.rm = TRUE)))
  
  fuera_rango <- sum(
    x < rango_min[i] | x > rango_max[i],
    na.rm = TRUE
  )
  cat("Fuera del rango orientativo [", rango_min[i], ",",
      rango_max[i], "]:", fuera_rango, "\n")
  
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  ric <- q[2] - q[1]
  atipicos <- x < q[1] - 1.5 * ric | x > q[2] + 1.5 * ric
  
  cat("Atípicos de Tukey:", sum(atipicos, na.rm = TRUE),
      sprintf("(%.2f%%)\n", 100 * mean(atipicos, na.rm = TRUE)))
  
  if (v != "precipitacion_mm") {
    r <- cor(x, dataset$precipitacion_mm, use = "complete.obs")
    cat("Pearson con precipitación:", round(r, 3), "\n")
  }
}

################################################################
# 2) Matriz de correlación en consola
################################################################

matriz_cor <- cor(dataset[, variables], use = "pairwise.complete.obs")
print(round(matriz_cor, 3))

################################################################
# 3) Panel de precipitación
################################################################

fila_max <- dataset[which.max(dataset$precipitacion_mm), ]
anio_max <- fila_max$anio
semana_max <- fila_max$semana

sub_max <- dataset[
  dataset$anio == anio_max & dataset$semana == semana_max,
  c("x", "y", "precipitacion_mm")
]

r_max <- rast(sub_max, type = "xyz", crs = "EPSG:4326")

promedio_celda <- aggregate(
  precipitacion_mm ~ id_celda + x + y,
  data = dataset, FUN = mean, na.rm = TRUE
)

r_precip <- rast(
  promedio_celda[, c("x", "y", "precipitacion_mm")],
  type = "xyz", crs = "EPSG:4326"
)

par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2),
    oma = c(0, 0, 2, 0), cex.main = 0.95)

hist(dataset$precipitacion_mm, breaks = 40, col = colores[1],
     border = "white", main = "Histograma",
     xlab = unidades[1], ylab = "Frecuencia")

boxplot(dataset$precipitacion_mm, col = colores[1],
        main = "Boxplot", ylab = unidades[1])

# Márgenes explícitos para evitar el recorte del título.
plot(r_precip,
     main = "Distribución espacial\npromedio",
     mar = c(4.5, 4.5, 4.5, 5),
     cex.main = 0.95,
     xlab = "Longitud", ylab = "Latitud",
     col = hcl.colors(50, "Blues", rev = TRUE))

# Se usa "banda" porque semana es la etiqueta del dataset,
# cuya correspondencia con la semana calendario requiere revisión.
plot(r_max,
     main = paste0("Evento máximo\nBanda ", semana_max, " | ", anio_max),
     mar = c(4.5, 4.5, 4.5, 5),
     cex.main = 0.95,
     xlab = "Longitud", ylab = "Latitud",
     col = hcl.colors(50, "Blues", rev = TRUE))

points(fila_max$x, fila_max$y, pch = 4, col = "red3", lwd = 2)

legend("bottomleft",
       legend = paste(round(fila_max$precipitacion_mm, 2), "mm"),
       pch = 4, col = "red3", bty = "n", cex = 0.8)

mtext("Precipitación | mm por periodo", outer = TRUE,
      side = 3, line = 0.4, font = 2)

cat("\nRegistro del máximo de precipitación:\n")
print(
  fila_max[, c("id_celda", "anio", "semana", "x", "y",
               "precipitacion_mm")],
  row.names = FALSE
)

################################################################
# 4) Panel individual de los seis predictores
################################################################

for (i in 2:length(variables)) {
  v <- variables[i]
  x <- dataset[[v]]
  
  # Relación resumida mediante 20 intervalos del predictor.
  ok <- is.finite(x) & is.finite(dataset$precipitacion_mm)
  bins <- cut(x[ok], breaks = 20, include.lowest = TRUE)
  
  promedio_por_bin <- tapply(dataset$precipitacion_mm[ok], bins, mean)
  centro_bin <- tapply(x[ok], bins, mean)
  bins_validos <- is.finite(centro_bin) & is.finite(promedio_por_bin)
  r <- cor(x[ok], dataset$precipitacion_mm[ok])
  
  datos_mapa <- data.frame(
    id_celda = dataset$id_celda,
    x = dataset$x, y = dataset$y, valor = x
  )
  
  promedio_celda <- aggregate(
    valor ~ id_celda + x + y,
    data = datos_mapa, FUN = mean, na.rm = TRUE
  )
  
  r_v <- rast(
    promedio_celda[, c("x", "y", "valor")],
    type = "xyz", crs = "EPSG:4326"
  )
  
  # Título del mapa según el tipo de variable.
  titulo_mapa <- if (v == "altitud_m") {
    "Altitud del terreno"
  } else if (grepl("^clim_", v)) {
    "Promedio climatológico\npor celda"
  } else {
    "Distribución espacial\npromedio"
  }
  
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2),
      oma = c(0, 0, 2, 0), cex.main = 0.95)
  
  hist(x, breaks = 40, col = colores[i], border = "white",
       main = "Histograma", xlab = unidades[i], ylab = "Frecuencia")
  
  boxplot(x, col = colores[i],
          main = "Boxplot", ylab = unidades[i])
  
  plot(centro_bin[bins_validos], promedio_por_bin[bins_validos],
       type = "o", pch = 16, cex = 0.6, col = colores[i], lwd = 2,
       xlab = unidades[i], ylab = "Precipitación media (mm/periodo)",
       main = paste0("Medias por intervalos\nPearson r = ", round(r, 3)))
  grid()
  
  # Márgenes explícitos para los mapas de todos los predictores.
  plot(r_v,
       main = titulo_mapa,
       mar = c(4.5, 4.5, 4.5, 5),
       cex.main = 0.95,
       xlab = "Longitud", ylab = "Latitud",
       col = hcl.colors(50, "Blues", rev = TRUE))
  
  mtext(paste(nombres[i], "|", unidades[i]), outer = TRUE,
        side = 3, line = 0.4, font = 2)
}

################################################################
# 5) Matriz de correlación como mapa de calor
################################################################

par(mfrow = c(1, 1), mar = c(8, 9, 4, 2),
    oma = c(0, 0, 0, 0), cex.main = 1)

n <- ncol(matriz_cor)
paleta <- colorRampPalette(c("#B2182B", "white", "#2166AC"))(201)

# Invertir filas para mostrar la matriz como una tabla.
image(1:n, 1:n, t(matriz_cor[n:1, ]), zlim = c(-1, 1),
      col = paleta, axes = FALSE, xlab = "", ylab = "",
      main = "Matriz de correlación\nPearson",
      xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5),
      xaxs = "i", yaxs = "i")

axis(1, at = 1:n, labels = nombres_matriz, las = 2,
     tick = FALSE, cex.axis = 0.85)

axis(2, at = 1:n, labels = rev(nombres_matriz), las = 2,
     tick = FALSE, cex.axis = 0.85)

abline(v = seq(0.5, n + 0.5), h = seq(0.5, n + 0.5),
       col = "grey85", lwd = 0.7)

for (fila in 1:n) {
  for (columna in 1:n) {
    valor <- matriz_cor[fila, columna]
    
    color_texto <- if (is.finite(valor) && abs(valor) >= 0.6) {
      "white"
    } else {
      "black"
    }
    
    text(columna, n - fila + 1,
         labels = if (is.finite(valor)) sprintf("%.2f", valor) else "NA",
         col = color_texto, cex = 0.9, font = 2)
  }
}

box()

################################################################
# Restablecer configuración gráfica
################################################################

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1),
    oma = c(0, 0, 0, 0), cex.main = 1.2)

cat("\nEDA terminado. Consulta los paneles anteriores con las flechas de Plots.\n")