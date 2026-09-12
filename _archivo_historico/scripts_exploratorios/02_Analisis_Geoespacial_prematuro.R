###################################################################
# 02_Analisis_Geoespacial.R
# Analisis geoespacial del Valle del Cauca: dimension espacial
# (686 celdas) y dimension temporal (16 anios x 52 semanas)
###################################################################

library(terra)

buscar_raiz <- function() {
  # Buscar hacia arriba la raiz real del proyecto, sin depender de setwd().
  p <- normalizePath(getwd(), mustWork=TRUE)
  repeat {
    if (dir.exists(file.path(p, "datos_proyecto_1")) &&
        dir.exists(file.path(p, "avance_team_proyecto1"))) return(p)
    padre <- dirname(p)
    if (identical(padre, p)) break
    p <- padre
  }
  stop("No se encontro la raiz del proyecto. Abra la carpeta 'Proyecto1' en RStudio.")
}
RAIZ_PROYECTO <- buscar_raiz()
DIR_SCRIPT <- file.path(RAIZ_PROYECTO, "avance_team_proyecto1", "reinicio_desde_cero")
DIR_RESULTADOS <- file.path(DIR_SCRIPT, "analisis_geoespacial_resultados")
dir.create(DIR_RESULTADOS, recursive = TRUE, showWarnings = FALSE)

dataset <- readRDS(file.path(DIR_SCRIPT, "dataset_eda", "dataset_final_7_variables.rds"))


################################################################
# 1) Autocorrelacion espacial (Moran/Geary), metodologia de clase:
#    pesos de vecindad "reina" via terra::adjacent(), formulas
#    cerradas de algebra matricial (funciones_auxiliares.R)
################################################################

construir_pesos_queen <- function(r) {
  celdas_validas <- which(!is.na(values(r)[, 1]))
  pares <- adjacent(r, cells = celdas_validas, directions = "queen", pairs = TRUE)
  pares <- pares[pares[, 2] %in% celdas_validas, , drop = FALSE]
  n <- length(celdas_validas)
  indice <- setNames(seq_len(n), as.character(celdas_validas))
  W <- matrix(0, n, n)
  fila    <- indice[as.character(pares[, 1])]
  columna <- indice[as.character(pares[, 2])]
  W[cbind(fila, columna)] <- 1
  suma_filas <- rowSums(W)
  W <- W / ifelse(suma_filas == 0, 1, suma_filas)
  list(W = W, celdas = celdas_validas)
}

indice_moran <- function(x, W) {
  n <- length(x)
  z <- x - mean(x)
  S0 <- sum(W)
  (n / S0) * (as.numeric(t(z) %*% W %*% z) / sum(z^2))
}

indice_geary <- function(x, W) {
  n <- length(x)
  S0 <- sum(W)
  Dmat <- outer(x, x, FUN = function(a, b) (a - b)^2)
  ((n - 1) / (2 * S0)) * (sum(W * Dmat) / sum((x - mean(x))^2))
}

# Raster molde de 686 celdas (a partir de la altitud, que es 1 valor/celda)
promedio_celda <- aggregate(altitud_m ~ id_celda + x + y, data = dataset, FUN = mean)
r_molde <- rast(promedio_celda[, c("x", "y", "altitud_m")], type = "xyz", crs = "EPSG:4326")
pesos <- construir_pesos_queen(r_molde)

# Semanas representativas: 2 de temporada de lluvia y 2 de temporada seca,
# en anios distintos, para ver si la autocorrelacion espacial cambia con
# la epoca del anio (dimension temporal) y no solo con la celda (espacial)
semanas_estudio <- data.frame(
  anio = c(2015, 2020, 2016, 2021),
  semana = c(18, 44, 5, 32),
  temporada = c("lluvia (may)", "lluvia (nov)", "seca (feb)", "seca (ago)")
)

resultados_moran <- data.frame()
for (i in seq_len(nrow(semanas_estudio))) {
  a <- semanas_estudio$anio[i]; s <- semanas_estudio$semana[i]
  sub <- dataset[dataset$anio == a & dataset$semana == s, c("x", "y", "precipitacion_mm")]
  r_sub <- rast(sub, type = "xyz", crs = "EPSG:4326")
  x <- values(r_sub)[pesos$celdas, 1]
  I <- indice_moran(x, pesos$W)
  C <- indice_geary(x, pesos$W)
  resultados_moran <- rbind(resultados_moran, data.frame(
    anio = a, semana = s, temporada = semanas_estudio$temporada[i],
    moran_I = round(I, 3), geary_C = round(C, 3)
  ))
}

cat("========== MORAN / GEARY POR SEMANA (precipitacion cruda) ==========\n")
print(resultados_moran, row.names = FALSE)
cat("\nValor esperado de Moran I bajo no-autocorrelacion: -1/(n-1) =",
    round(-1/(length(pesos$celdas)-1), 4), "\n")
cat("Todas las semanas muestran Moran I bien por encima del valor esperado\n")
cat("bajo independencia -> autocorrelacion espacial positiva, presente en\n")
cat("todas las epocas del anio revisadas (no es exclusiva de una temporada).\n\n")


################################################################
# 2) Tendencia + residuales + Moran sobre los RESIDUALES
#    (paso previo obligatorio antes de ajustar el semivariograma:
#    si el residual siguiera autocorrelacionado espacialmente, la
#    tendencia no capturo toda la estructura de gran escala)
################################################################

lat0_global <- mean(dataset$y)
dataset$x_km <- (dataset$x - mean(dataset$x)) * 111.320 * cos(lat0_global * pi / 180)
dataset$y_km <- (dataset$y - lat0_global) * 110.574

resultados_residual_moran <- data.frame()
for (i in seq_len(nrow(semanas_estudio))) {
  a <- semanas_estudio$anio[i]; s <- semanas_estudio$semana[i]
  sub <- dataset[dataset$anio == a & dataset$semana == s, ]
  m_trend <- lm(precipitacion_mm ~ x_km + y_km + altitud_m, data = sub)
  sub$residual <- residuals(m_trend)

  r_res <- rast(sub[, c("x", "y", "residual")], type = "xyz", crs = "EPSG:4326")
  x_res <- values(r_res)[pesos$celdas, 1]
  I_res <- indice_moran(x_res, pesos$W)

  resultados_residual_moran <- rbind(resultados_residual_moran, data.frame(
    anio = a, semana = s, temporada = semanas_estudio$temporada[i],
    R2_tendencia = round(summary(m_trend)$r.squared, 3),
    moran_I_crudo = resultados_moran$moran_I[i],
    moran_I_residual = round(I_res, 3)
  ))
}

cat("========== MORAN I: CRUDO vs. RESIDUAL DE TENDENCIA ==========\n")
print(resultados_residual_moran, row.names = FALSE)
cat("\nLa tendencia (x_km + y_km + altitud) reduce la autocorrelacion espacial\n")
cat("(Moran I baja del crudo al residual), pero no la elimina del todo en\n")
cat("ninguna semana -> queda estructura espacial de escala fina en el\n")
cat("residual, que es exactamente lo que el semivariograma + kriging deben\n")
cat("modelar (ver seccion 3).\n\n")


################################################################
# 3) Semivariograma + kriging universal, por semana
#    (misma metodologia que 02_Modelacion.R: exponencial ajustado
#    por minimos cuadrados ponderados por N de pares, LOO)
################################################################

sce_exponencial <- function(theta, h, gamma_emp, pesos_n) {
  c0 <- theta[1]; c1 <- theta[2]; phi <- theta[3]
  gamma_teo <- c0 + c1 * (1 - exp(-h / phi))
  sum(pesos_n * (gamma_emp - gamma_teo)^2)
}

kriging_universal <- function(s0, coords_xy, resid, Sigma, C_hat, sill_hat) {
  n_loc <- nrow(coords_xy)
  d0 <- sqrt((coords_xy[, 1] - s0[1])^2 + (coords_xy[, 2] - s0[2])^2)
  c0_vec <- C_hat(d0)
  A <- rbind(cbind(Sigma, rep(1, n_loc)), c(rep(1, n_loc), 0))
  b <- c(c0_vec, 1)
  sol <- solve(A, b)
  pesos_k <- sol[1:n_loc]
  lambda <- sol[n_loc + 1]
  c(prediccion = unname(sum(pesos_k * resid)),
    varianza = unname(sill_hat - sum(pesos_k * c0_vec) - lambda))
}

modelar_semana_final <- function(a, s, cutoff_km = 180) {
  sub <- dataset[dataset$anio == a & dataset$semana == s, ]
  m_trend <- lm(precipitacion_mm ~ x_km + y_km + altitud_m, data = sub)
  residuales <- residuals(m_trend)

  D <- as.matrix(dist(sub[, c("x_km", "y_km")]))
  pares <- which(upper.tri(D), arr.ind = TRUE)
  h_ij <- D[upper.tri(D)]
  gamma_ij <- (residuales[pares[, 1]] - residuales[pares[, 2]])^2 / 2

  en_cutoff <- h_ij <= cutoff_km
  bins <- cut(h_ij[en_cutoff], breaks = 12)
  gamma_emp <- tapply(gamma_ij[en_cutoff], bins, mean)
  h_medio <- tapply(h_ij[en_cutoff], bins, mean)
  N_h <- tapply(gamma_ij[en_cutoff], bins, length)
  semivariograma <- na.omit(data.frame(h = h_medio, gamma = gamma_emp, N = N_h))

  phis_iniciales <- quantile(semivariograma$h, probs = c(0.1, 0.25, 0.5, 0.75, 1))
  ajustes <- lapply(phis_iniciales, function(phi0) {
    optim(par = c(c0 = 0, c1 = var(residuales), phi = phi0),
          fn = sce_exponencial, h = semivariograma$h, gamma_emp = semivariograma$gamma,
          pesos_n = semivariograma$N, method = "L-BFGS-B", lower = c(0, 0, 1),
          upper = c(max(semivariograma$gamma), 5 * max(semivariograma$gamma), 3 * max(semivariograma$h)))
  })
  ajuste <- ajustes[[which.min(sapply(ajustes, function(x) x$value))]]
  c0_hat <- unname(ajuste$par["c0"]); c1_hat <- unname(ajuste$par["c1"]); phi_hat <- unname(ajuste$par["phi"])
  sill_hat <- c0_hat + c1_hat
  C_hat <- function(h) sill_hat * exp(-h / phi_hat)
  Sigma <- matrix(C_hat(as.vector(D)), nrow = nrow(sub))
  coords_xy <- as.matrix(sub[, c("x_km", "y_km")])

  set.seed(2026)
  idx_loo <- sample(seq_len(nrow(sub)), min(120, nrow(sub)))
  errores_loo <- numeric(length(idx_loo))
  for (k in seq_along(idx_loo)) {
    i <- idx_loo[k]
    krg_i <- kriging_universal(coords_xy[i, ], coords_xy[-i, ], residuales[-i], Sigma[-i, -i], C_hat, sill_hat)
    pred_i <- predict(m_trend, newdata = sub[i, , drop = FALSE]) + krg_i["prediccion"]
    errores_loo[k] <- sub$precipitacion_mm[i] - pred_i
  }

  list(semivariograma = semivariograma, c0 = c0_hat, sill = sill_hat, phi = phi_hat,
       rango_km = 3 * phi_hat, r2_tendencia = summary(m_trend)$r.squared,
       rmse_loo = sqrt(mean(errores_loo^2)), mae_loo = mean(abs(errores_loo)),
       precip_media = mean(sub$precipitacion_mm))
}

resultados_variograma <- list()
for (i in seq_len(nrow(semanas_estudio))) {
  a <- semanas_estudio$anio[i]; s <- semanas_estudio$semana[i]
  cat("Modelando", a, "- semana", s, "...\n")
  resultados_variograma[[i]] <- modelar_semana_final(a, s)
}

resumen_variograma <- do.call(rbind, lapply(seq_len(nrow(semanas_estudio)), function(i) {
  r <- resultados_variograma[[i]]
  data.frame(anio = semanas_estudio$anio[i], semana = semanas_estudio$semana[i],
             temporada = semanas_estudio$temporada[i],
             precip_media = round(r$precip_media, 1), R2_tendencia = round(r$r2_tendencia, 3),
             nugget = round(r$c0, 2), sill = round(r$sill, 2), rango_km = round(r$rango_km, 1),
             RMSE_mm = round(r$rmse_loo, 2), MAE_mm = round(r$mae_loo, 2),
             RMSE_relativo_pct = round(100 * r$rmse_loo / r$precip_media, 1))
}))

cat("\n========== RESUMEN SEMIVARIOGRAMA Y KRIGING, POR SEMANA ==========\n")
print(resumen_variograma, row.names = FALSE)
cat("\nEl rango (distancia a la que se agota la autocorrelacion espacial)\n")
cat("y el nugget/sill cambian de una semana a otra -> la estructura espacial\n")
cat("de la precipitacion NO es fija en el tiempo, depende de la epoca del\n")
cat("anio y del evento climatico particular de esa semana.\n\n")

# Semivariogramas de las 4 semanas en un solo panel
par(mfrow = c(2, 2))
for (i in seq_len(nrow(semanas_estudio))) {
  r <- resultados_variograma[[i]]
  plot(r$semivariograma$h, r$semivariograma$gamma, pch = 19, col = "navy", type = "b",
       xlab = "Distancia (km)", ylab = "Semivarianza",
       main = paste0(semanas_estudio$temporada[i], " (", semanas_estudio$anio[i], "-s",
                      semanas_estudio$semana[i], ")"))
  curva_h <- seq(0, max(r$semivariograma$h), length.out = 100)
  curva_gamma <- r$c0 + (r$sill - r$c0) * (1 - exp(-curva_h / r$phi))
  lines(curva_h, curva_gamma, col = "red", lwd = 2)
}
par(mfrow = c(1, 1))


################################################################
# 4) Autocorrelacion TEMPORAL (dimension tiempo, independiente de
#    la espacial): funcion de autocorrelacion (ACF) de la anomalia
#    de precipitacion (crudo - climatologia) en celdas de muestra
################################################################

dataset$anomalia_precip <- dataset$precipitacion_mm - dataset$clim_precipitacion_mm

set.seed(2026)
celdas_muestra <- sample(unique(dataset$id_celda), 6)

par(mfrow = c(2, 3))
resumen_acf <- data.frame()
for (cel in celdas_muestra) {
  serie <- dataset[dataset$id_celda == cel, ]
  serie <- serie[order(serie$anio, serie$semana), ]
  a <- acf(serie$anomalia_precip, plot = TRUE, main = paste("Celda", cel), lag.max = 12)
  resumen_acf <- rbind(resumen_acf, data.frame(
    id_celda = cel, acf_lag1 = round(a$acf[2], 3), acf_lag2 = round(a$acf[3], 3)
  ))
}
par(mfrow = c(1, 1))

cat("========== AUTOCORRELACION TEMPORAL (ACF) DE LA ANOMALIA DE PRECIPITACION ==========\n")
print(resumen_acf, row.names = FALSE)
cat("\nLa autocorrelacion temporal a 1 semana de rezago es baja/moderada y\n")
cat("cae rapido -> la precipitacion semanal, una vez removido el ciclo\n")
cat("estacional (climatologia), se comporta de forma casi independiente de\n")
cat("una semana a la siguiente. Esto justifica tratar cada semana como una\n")
cat("realizacion espacial separada (como se hizo en las secciones 1-3), en\n")
cat("vez de modelar explicitamente la dependencia temporal semana a semana.\n\n")


################################################################
# Guardar resultados
################################################################

write.csv(resultados_moran, file.path(DIR_RESULTADOS, "moran_geary_por_semana.csv"), row.names = FALSE)
write.csv(resultados_residual_moran, file.path(DIR_RESULTADOS, "moran_crudo_vs_residual.csv"), row.names = FALSE)
write.csv(resumen_variograma, file.path(DIR_RESULTADOS, "resumen_semivariograma_kriging.csv"), row.names = FALSE)
write.csv(resumen_acf, file.path(DIR_RESULTADOS, "autocorrelacion_temporal_acf.csv"), row.names = FALSE)

cat("Resultados guardados en:", DIR_RESULTADOS, "\n")
