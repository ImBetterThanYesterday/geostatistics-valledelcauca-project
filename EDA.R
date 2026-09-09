###################################################################
# EDA por variable - dataset final (7 variables)
# Proyecto 1 - Analitica de Datos
###################################################################

dataset <- readRDS("dataset_eda/dataset_final_7_variables.rds")
str(dataset)


################################################################
# 0) Verificacion de cobertura temporal: meses por anio
################################################################

dataset$mes_aprox <- pmin(12, ceiling(dataset$semana * 7 / 30.4))
tabla_meses <- table(dataset$anio, dataset$mes_aprox) > 0
tabla_meses

apply(tabla_meses, 1, sum)

# Los 16 anios (2010-2025) tienen los 12 meses cubiertos -> no hay
# anios incompletos, el soporte temporal esta completo.


################################################################
# 1) Estadistica descriptiva, nulos, atipicos e incoherentes
#    por variable
################################################################

variables  <- c("precipitacion_mm", "temperatura_bilineal_c", "radiacion_bilineal_mj_m2_dia",
                "altitud_m", "clim_precipitacion_mm", "clim_temperatura_bilineal_c",
                "clim_radiacion_bilineal_mj_m2_dia")
rango_min <- c(0, -10, 0, -50, 0, -10, 0)
rango_max <- c(700, 45, 45, 6000, 700, 45, 45)
colores   <- c("steelblue", "firebrick", "orange", "darkgreen",
               "steelblue", "firebrick", "orange")

for (i in seq_along(variables)) {

  v <- variables[i]
  x <- dataset[[v]]

  cat("\n============================================================\n")
  cat(v, "\n")
  cat("============================================================\n")

  # Estadistica descriptiva (incluye min y max)
  print(summary(x))
  cat("sd:", sd(x), "\n")

  # Nulos
  cat("Nulos:", sum(is.na(x)), "\n")

  # Incoherentes (fuera de rango fisico)
  incoherentes <- sum(x < rango_min[i] | x > rango_max[i], na.rm = TRUE)
  cat("Incoherentes (fuera de [", rango_min[i], ",", rango_max[i], "]):", incoherentes, "\n")

  # Atipicos (regla de Tukey)
  q <- quantile(x, c(0.25, 0.75))
  ric <- q[2] - q[1]
  atipicos <- x < (q[1] - 1.5 * ric) | x > (q[2] + 1.5 * ric)
  cat("Atipicos (Tukey):", sum(atipicos), sprintf("(%.2f%%)\n", 100 * mean(atipicos)))

  # Histograma
  hist(x, breaks = 40, col = colores[i], main = paste("Histograma -", v), xlab = v)

  # Boxplot (caja de bigotes)
  boxplot(x, main = paste("Boxplot -", v), ylab = v)

  # Dispersion contra precipitacion (no se hace para precipitacion
  # misma, seria la variable contra si misma)
  if (v != "precipitacion_mm") {
    plot(x, dataset$precipitacion_mm, pch = 20, col = colores[i],
         xlab = v, ylab = "precipitacion_mm",
         main = paste("Precipitacion vs.", v))
    cat("Correlacion con precipitacion:", round(cor(x, dataset$precipitacion_mm), 3), "\n")
  }
}


################################################################
# 2) Imagen combinada: las 6 variables restantes vs. precipitacion
################################################################

par(mfrow = c(2, 3))
plot(dataset$temperatura_bilineal_c, dataset$precipitacion_mm, pch = 20, col = "firebrick",
     xlab = "temperatura_bilineal_c", ylab = "precipitacion_mm", main = "Precip. vs. Temperatura")
plot(dataset$radiacion_bilineal_mj_m2_dia, dataset$precipitacion_mm, pch = 20, col = "orange",
     xlab = "radiacion_bilineal_mj_m2_dia", ylab = "precipitacion_mm", main = "Precip. vs. Radiacion")
plot(dataset$altitud_m, dataset$precipitacion_mm, pch = 20, col = "darkgreen",
     xlab = "altitud_m", ylab = "precipitacion_mm", main = "Precip. vs. Altitud")
plot(dataset$clim_precipitacion_mm, dataset$precipitacion_mm, pch = 20, col = "steelblue",
     xlab = "clim_precipitacion_mm", ylab = "precipitacion_mm", main = "Precip. vs. Clim. precip.")
plot(dataset$clim_temperatura_bilineal_c, dataset$precipitacion_mm, pch = 20, col = "firebrick",
     xlab = "clim_temperatura_bilineal_c", ylab = "precipitacion_mm", main = "Precip. vs. Clim. temp.")
plot(dataset$clim_radiacion_bilineal_mj_m2_dia, dataset$precipitacion_mm, pch = 20, col = "orange",
     xlab = "clim_radiacion_bilineal_mj_m2_dia", ylab = "precipitacion_mm", main = "Precip. vs. Clim. rad.")
par(mfrow = c(1, 1))

# La relacion mas fuerte a simple vista es con la climatologia de
# precipitacion (0.71), seguida de climatologia de temperatura (0.50)
# y temperatura cruda (0.45). Altitud tiene una relacion moderada y
# negativa (-0.39). Radiacion (cruda y climatologica) es la mas
# debil de las 6, con nubes de puntos bastante dispersas.


################################################################
# 3) Matriz de correlacion
################################################################

matriz_cor <- cor(dataset[, variables])
round(matriz_cor, 2)

# Precipitacion correlaciona mas fuerte con su propia climatologia
# (0.71) y con temperatura (0.45); negativo con altitud (-0.39) y
# radiacion (-0.22). Temperatura y altitud tienen la correlacion mas
# fuerte de toda la matriz (-0.66), el gradiente termico.
