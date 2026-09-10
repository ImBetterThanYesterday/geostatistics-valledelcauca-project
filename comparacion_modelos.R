# Proyecto 1 - Analitica de Datos
# Fase 2: Comparacion de escenarios de modelacion espacial
# Validacion cruzada por bloques espaciales

library(terra)
library(ggplot2)
library(gridExtra)

ruta_dataset <- "C:/Users/Administrador/Documents/01_PROYECTO_ANALITICA/dataset_eda/dataset_final_7_variables.rds"

set.seed(2026)

lado_bloque_km <- 35
n_folds <- 5
n_bins <- 15
prop_cutoff <- 0.35


# 1. Carga y preparacion

datos <- readRDS(ruta_dataset)

datos$x_km <- datos$x * 111.32 * cos(datos$y * pi / 180)
datos$y_km <- datos$y * 110.57

datos$temperatura_c <- datos$temperatura_bilineal_c
datos$radiacion_mj <- datos$radiacion_bilineal_mj_m2_dia

semanas_estudio <- data.frame(
  anio   = c(2011, 2012, 2014, 2015, 2016, 2018, 2020, 2021, 2023, 2024),
  semana = c(8,    18,   43,   26,   4,    41,   15,   45,   9,    10)
)

vars_modelo <- c("x_km", "y_km", "altitud_m", "precipitacion_mm",
                 "temperatura_c", "radiacion_mj")


# 2. Funciones geoestadisticas

semivariograma_empirico <- function(coords, valores, n_bins, prop_cutoff) {
  d <- as.matrix(dist(coords))
  idx <- upper.tri(d)
  h <- d[idx]
  
  dif <- outer(valores, valores, "-")
  semiv <- (dif[idx]^2) / 2
  
  cutoff <- max(h) * prop_cutoff
  sel <- h <= cutoff & h > 0
  h <- h[sel]
  semiv <- semiv[sel]
  
  bins <- cut(h, breaks = seq(0, cutoff, length.out = n_bins + 1), include.lowest = TRUE)
  
  vg <- data.frame(
    dist  = as.numeric(tapply(h, bins, mean)),
    gamma = as.numeric(tapply(semiv, bins, mean)),
    n     = as.numeric(tapply(h, bins, length))
  )
  vg[!is.na(vg$dist), ]
}

modelo_exponencial <- function(h, nugget, meseta, rango) {
  nugget + meseta * (1 - exp(-h / rango))
}

ajustar_variograma <- function(vg) {
  sce <- function(p) {
    pred <- modelo_exponencial(vg$dist, abs(p[1]), abs(p[2]), abs(p[3]))
    sum(vg$n * (vg$gamma - pred)^2)
  }
  
  ini <- c(min(vg$gamma), max(vg$gamma) - min(vg$gamma), max(vg$dist) / 3)
  fit <- optim(ini, sce, method = "Nelder-Mead",
               control = list(maxit = 5000, reltol = 1e-10))
  
  list(nugget = abs(fit$par[1]), meseta = abs(fit$par[2]), rango = abs(fit$par[3]))
}

kriging_ordinario <- function(coords_ent, valores_ent, coords_pred, par) {
  n <- nrow(coords_ent)
  
  G <- modelo_exponencial(as.matrix(dist(coords_ent)), par$nugget, par$meseta, par$rango)
  diag(G) <- 0
  
  A <- matrix(1, n + 1, n + 1)
  A[1:n, 1:n] <- G
  A[n + 1, n + 1] <- 0
  
  A_inv <- tryCatch(solve(A), error = function(e) MASS::ginv(A))
  
  pred <- numeric(nrow(coords_pred))
  for (i in seq_len(nrow(coords_pred))) {
    d_pe <- sqrt((coords_ent[, 1] - coords_pred[i, 1])^2 +
                   (coords_ent[, 2] - coords_pred[i, 2])^2)
    b <- c(modelo_exponencial(d_pe, par$nugget, par$meseta, par$rango), 1)
    lambda <- A_inv %*% b
    pred[i] <- sum(lambda[1:n] * valores_ent)
  }
  pred
}

metricas <- function(obs, pred) {
  e <- obs - pred
  c(rmse = sqrt(mean(e^2)),
    mae  = mean(abs(e)),
    r2   = 1 - sum(e^2) / sum((obs - mean(obs))^2))
}


# 3. Escenarios

escenarios <- list(
  M0 = NULL,
  M1 = precipitacion_mm ~ x_km + y_km,
  M2 = precipitacion_mm ~ altitud_m,
  M3 = precipitacion_mm ~ x_km + y_km + altitud_m,
  M4 = precipitacion_mm ~ altitud_m + temperatura_c + radiacion_mj
)


# 4. Asignacion de bloques espaciales

asignar_bloques <- function(d, lado) {
  bx <- floor((d$x_km - min(d$x_km)) / lado)
  by <- floor((d$y_km - min(d$y_km)) / lado)
  paste(bx, by, sep = "_")
}

asignar_folds <- function(bloques, k) {
  b <- unique(bloques)
  asignacion <- sample(rep(seq_len(k), length.out = length(b)))
  names(asignacion) <- b
  asignacion[bloques]
}


# 5. Validacion por semana

ajustar_semana <- function(anio_i, semana_i) {
  d <- datos[datos$anio == anio_i & datos$semana == semana_i, ]
  d <- d[complete.cases(d[, vars_modelo]), ]
  if (nrow(d) < 100) return(NULL)
  
  d$bloque <- asignar_bloques(d, lado_bloque_km)
  d$fold <- asignar_folds(d$bloque, n_folds)
  
  filas <- list()
  
  for (nombre in names(escenarios)) {
    f <- escenarios[[nombre]]
    
    obs <- c(); pred_k <- c(); pred_reg <- c()
    
    for (kf in seq_len(n_folds)) {
      ent <- d[d$fold != kf, ]
      pru <- d[d$fold == kf, ]
      if (nrow(pru) < 5 || nrow(ent) < 50) next
      
      coords_ent <- as.matrix(ent[, c("x_km", "y_km")])
      coords_pru <- as.matrix(pru[, c("x_km", "y_km")])
      
      if (is.null(f)) {
        residual_ent <- ent$precipitacion_mm
        tendencia <- rep(0, nrow(pru))
      } else {
        m <- lm(f, data = ent)
        residual_ent <- residuals(m)
        tendencia <- predict(m, newdata = pru)
      }
      
      vg <- semivariograma_empirico(coords_ent, residual_ent, n_bins, prop_cutoff)
      par <- ajustar_variograma(vg)
      
      correccion <- kriging_ordinario(coords_ent, residual_ent, coords_pru, par)
      
      obs <- c(obs, pru$precipitacion_mm)
      pred_k <- c(pred_k, tendencia + correccion)
      pred_reg <- c(pred_reg, if (is.null(f)) rep(mean(ent$precipitacion_mm), nrow(pru)) else tendencia)
    }
    
    m_k <- metricas(obs, pred_k)
    m_r <- metricas(obs, pred_reg)
    
    filas[[nombre]] <- data.frame(
      anio = anio_i, semana = semana_i, modelo = nombre,
      rmse = m_k["rmse"], mae = m_k["mae"], r2 = m_k["r2"],
      rmse_sin_kriging = m_r["rmse"],
      mae_sin_kriging = m_r["mae"],
      r2_sin_kriging = m_r["r2"],
      n_bloques = length(unique(d$bloque)),
      n_total = nrow(d)
    )
  }
  
  do.call(rbind, filas)
}


# 6. Ejecucion

resultados <- list()

for (k in seq_len(nrow(semanas_estudio))) {
  a <- semanas_estudio$anio[k]
  s <- semanas_estudio$semana[k]
  cat("Procesando", a, "semana", s, "\n")
  resultados[[k]] <- ajustar_semana(a, s)
}

resultados <- do.call(rbind, resultados)
rownames(resultados) <- NULL

cat("\nBloques espaciales por semana:", unique(resultados$n_bloques), "\n")


# 7. Resumen comparativo

resumen <- data.frame(
  modelo = names(escenarios),
  rmse_medio = tapply(resultados$rmse, resultados$modelo, mean)[names(escenarios)],
  mae_medio  = tapply(resultados$mae, resultados$modelo, mean)[names(escenarios)],
  r2_medio   = tapply(resultados$r2, resultados$modelo, mean)[names(escenarios)],
  rmse_sin_kriging = tapply(resultados$rmse_sin_kriging, resultados$modelo, mean)[names(escenarios)]
)
resumen$mejora_kriging_pct <- 100 * (resumen$rmse_sin_kriging - resumen$rmse_medio) /
  resumen$rmse_sin_kriging

cat("\nResumen comparativo con validacion por bloques\n\n")
print(round(resumen[, -1], 3))

mejor <- resumen$modelo[which.min(resumen$rmse_medio)]
cat("\nMejor escenario:", mejor, "\n")


# 8. Comparacion directa contra M0

base_m0 <- resultados[resultados$modelo == "M0", c("anio", "semana", "rmse")]
names(base_m0)[3] <- "rmse_m0"
comp <- merge(resultados, base_m0, by = c("anio", "semana"))
comp$dif_vs_m0 <- comp$rmse - comp$rmse_m0

cat("\nDiferencia media de RMSE frente a M0 (negativo = mejor que M0)\n")
print(round(tapply(comp$dif_vs_m0, comp$modelo, mean)[names(escenarios)], 3))

cat("\nSemanas en que cada modelo supera a M0 (de", nrow(semanas_estudio), ")\n")
print(tapply(comp$dif_vs_m0 < 0, comp$modelo, sum)[names(escenarios)])


# 9. Graficos

g1 <- ggplot(resultados, aes(x = modelo, y = rmse, fill = modelo)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "RMSE por escenario (bloques espaciales)", x = NULL, y = "RMSE (mm)") +
  theme_minimal()

g2 <- ggplot(resultados, aes(x = modelo, y = mae, fill = modelo)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "MAE por escenario", x = NULL, y = "MAE (mm)") +
  theme_minimal()

grid.arrange(g1, g2, ncol = 2)

resultados$periodo <- paste0(resultados$anio, "_s", resultados$semana)

g3 <- ggplot(resultados, aes(x = periodo, y = rmse, group = modelo, color = modelo)) +
  geom_line() + geom_point() +
  labs(title = "RMSE por semana y escenario", x = NULL, y = "RMSE (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(g3)

comp$periodo <- paste0(comp$anio, "_s", comp$semana)
g4 <- ggplot(comp[comp$modelo != "M0", ],
             aes(x = periodo, y = dif_vs_m0, fill = modelo)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0) +
  labs(title = "Diferencia de RMSE frente a M0", x = NULL, y = "RMSE - RMSE(M0)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(g4)


# 10. Variograma del mejor escenario

d_ref <- datos[datos$anio == semanas_estudio$anio[1] &
                 datos$semana == semanas_estudio$semana[1], ]
d_ref <- d_ref[complete.cases(d_ref[, vars_modelo]), ]

f_mejor <- escenarios[[mejor]]
res_ref <- if (is.null(f_mejor)) d_ref$precipitacion_mm else residuals(lm(f_mejor, data = d_ref))

vg_ref <- semivariograma_empirico(as.matrix(d_ref[, c("x_km", "y_km")]), res_ref, n_bins, prop_cutoff)
par_ref <- ajustar_variograma(vg_ref)

curva <- data.frame(dist = seq(0, max(vg_ref$dist), length.out = 200))
curva$gamma <- modelo_exponencial(curva$dist, par_ref$nugget, par_ref$meseta, par_ref$rango)

g5 <- ggplot(vg_ref, aes(x = dist, y = gamma)) +
  geom_point(aes(size = n)) +
  geom_line(data = curva, color = "firebrick", linewidth = 1) +
  labs(title = paste("Semivariograma de residuales -", mejor),
       x = "Distancia (km)", y = "Semivarianza") +
  theme_minimal()

print(g5)

cat("\nParametros del variograma\n")
cat("nugget:", round(par_ref$nugget, 3),
    "| meseta:", round(par_ref$meseta, 3),
    "| rango:", round(par_ref$rango, 2), "km\n")


# 11. Mapa de bloques usado en la validacion

d_bloques <- datos[datos$anio == semanas_estudio$anio[1] &
                     datos$semana == semanas_estudio$semana[1], ]
d_bloques <- d_bloques[complete.cases(d_bloques[, vars_modelo]), ]
d_bloques$bloque <- asignar_bloques(d_bloques, lado_bloque_km)
d_bloques$fold <- factor(asignar_folds(d_bloques$bloque, n_folds))

g6 <- ggplot(d_bloques, aes(x = x, y = y, color = fold)) +
  geom_point(size = 1.2) +
  coord_equal() +
  labs(title = paste0("Particion en bloques espaciales (", lado_bloque_km, " km)"),
       x = "Longitud", y = "Latitud") +
  theme_minimal()

print(g6)




# =====================================================================
# Fase 3: Prediccion y validacion del modelo seleccionado (M0)
# =====================================================================

semana_pred <- list(anio = 2011, semana = 8)
resolucion_pred <- 0.02
radio_mascara_km <- 6


# 12. Kriging ordinario con varianza de prediccion

kriging_ordinario_var <- function(coords_ent, valores_ent, coords_pred, par) {
  n <- nrow(coords_ent)
  
  G <- modelo_exponencial(as.matrix(dist(coords_ent)), par$nugget, par$meseta, par$rango)
  diag(G) <- 0
  
  A <- matrix(1, n + 1, n + 1)
  A[1:n, 1:n] <- G
  A[n + 1, n + 1] <- 0
  A_inv <- solve(A)
  
  pred <- numeric(nrow(coords_pred))
  varianza <- numeric(nrow(coords_pred))
  
  for (i in seq_len(nrow(coords_pred))) {
    d_pe <- sqrt((coords_ent[, 1] - coords_pred[i, 1])^2 +
                   (coords_ent[, 2] - coords_pred[i, 2])^2)
    b <- c(modelo_exponencial(d_pe, par$nugget, par$meseta, par$rango), 1)
    lambda <- A_inv %*% b
    pred[i] <- sum(lambda[1:n] * valores_ent)
    varianza[i] <- sum(lambda * b)
  }
  
  list(pred = pred, var = varianza)
}


# 13. Validacion dejando uno fuera del modelo final

d_val <- datos[datos$anio == semana_pred$anio & datos$semana == semana_pred$semana, ]
d_val <- d_val[complete.cases(d_val[, vars_modelo]), ]

coords_val <- as.matrix(d_val[, c("x_km", "y_km")])
z_val <- d_val$precipitacion_mm

vg_val <- semivariograma_empirico(coords_val, z_val, n_bins, prop_cutoff)
par_val <- ajustar_variograma(vg_val)

cat("\nVariograma del modelo final\n")
cat("nugget:", round(par_val$nugget, 3),
    "| meseta:", round(par_val$meseta, 2),
    "| rango:", round(par_val$rango, 2), "km\n")

n_val <- nrow(d_val)
pred_loo <- numeric(n_val)
var_loo <- numeric(n_val)

for (i in seq_len(n_val)) {
  k <- kriging_ordinario_var(coords_val[-i, , drop = FALSE], z_val[-i],
                             coords_val[i, , drop = FALSE], par_val)
  pred_loo[i] <- k$pred
  var_loo[i] <- k$var
}

d_val$pred_loo <- pred_loo
d_val$error_loo <- z_val - pred_loo
d_val$var_loo <- var_loo
d_val$error_estandarizado <- d_val$error_loo / sqrt(pmax(var_loo, 1e-8))

m_loo <- metricas(z_val, pred_loo)

cat("\nValidacion dejando uno fuera (n =", n_val, ")\n")
cat("RMSE:", round(m_loo["rmse"], 3), "mm\n")
cat("MAE :", round(m_loo["mae"], 3), "mm\n")
cat("R2  :", round(m_loo["r2"], 3), "\n")
cat("Sesgo medio:", round(mean(d_val$error_loo), 3), "mm\n")
cat("Media del error estandarizado:", round(mean(d_val$error_estandarizado), 3),
    "(esperado 0)\n")
cat("Varianza del error estandarizado:", round(var(d_val$error_estandarizado), 3),
    "(esperado 1)\n")


# 14. Graficos de validacion

g_disp <- ggplot(d_val, aes(x = pred_loo, y = precipitacion_mm)) +
  geom_point(alpha = 0.5, size = 1.3) +
  geom_abline(intercept = 0, slope = 1, color = "firebrick", linewidth = 0.8) +
  coord_equal() +
  labs(title = "Observado vs predicho (LOO)",
       x = "Predicho (mm)", y = "Observado (mm)") +
  theme_minimal()

g_res <- ggplot(d_val, aes(x = pred_loo, y = error_loo)) +
  geom_point(alpha = 0.5, size = 1.3) +
  geom_hline(yintercept = 0, color = "firebrick") +
  labs(title = "Residuales de validacion",
       x = "Predicho (mm)", y = "Error (mm)") +
  theme_minimal()

g_hist <- ggplot(d_val, aes(x = error_loo)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Distribucion del error", x = "Error (mm)", y = "Frecuencia") +
  theme_minimal()

g_mapa_err <- ggplot(d_val, aes(x = x, y = y, color = error_loo)) +
  geom_point(size = 1.4) +
  scale_color_gradient2(low = "#b2182b", mid = "grey90", high = "#2166ac", midpoint = 0) +
  coord_equal() +
  labs(title = "Distribucion espacial del error", x = "Longitud", y = "Latitud",
       color = "Error (mm)") +
  theme_minimal()

grid.arrange(g_disp, g_res, g_hist, g_mapa_err, ncol = 2)


# 15. Moran del error de validacion

matriz_pesos <- function(coords, umbral_km) {
  d <- as.matrix(dist(coords))
  W <- ifelse(d > 0 & d <= umbral_km, 1, 0)
  W
}

indice_moran <- function(x, W) {
  n <- length(x)
  z <- x - mean(x)
  s0 <- sum(W)
  (n / s0) * (t(z) %*% W %*% z) / sum(z^2)
}

W_val <- matriz_pesos(coords_val, 8)
I_error <- indice_moran(d_val$error_loo, W_val)
I_obs <- indice_moran(z_val, W_val)

cat("\nIndice de Moran\n")
cat("Precipitacion observada:", round(I_obs, 4), "\n")
cat("Error de validacion    :", round(I_error, 4), "\n")
cat("Valor esperado bajo H0 :", round(-1 / (n_val - 1), 4), "\n")


# 16. Prediccion sobre grilla fina

rango_x <- range(d_val$x)
rango_y <- range(d_val$y)

grilla <- expand.grid(
  x = seq(rango_x[1], rango_x[2], by = resolucion_pred),
  y = seq(rango_y[1], rango_y[2], by = resolucion_pred)
)

d_cel <- as.matrix(dist(rbind(
  as.matrix(grilla),
  as.matrix(d_val[, c("x", "y")])
)))[seq_len(nrow(grilla)), -seq_len(nrow(grilla))]

dist_min_km <- apply(d_cel, 1, min) * 111
grilla <- grilla[dist_min_km <= radio_mascara_km, ]

grilla$x_km <- grilla$x * 111.32 * cos(grilla$y * pi / 180)
grilla$y_km <- grilla$y * 110.57

cat("\nPuntos de la grilla de prediccion:", nrow(grilla), "\n")

k_grilla <- kriging_ordinario_var(coords_val, z_val,
                                  as.matrix(grilla[, c("x_km", "y_km")]), par_val)

grilla$precipitacion <- k_grilla$pred
grilla$desv <- sqrt(pmax(k_grilla$var, 0))

cat("Rango predicho:", round(min(grilla$precipitacion), 2), "-",
    round(max(grilla$precipitacion), 2), "mm\n")
cat("Rango observado:", round(min(z_val), 2), "-", round(max(z_val), 2), "mm\n")


# 17. Mapas de prediccion

etiqueta <- paste0(semana_pred$anio, " semana ", semana_pred$semana)

g_obs <- ggplot(d_val, aes(x = x, y = y, color = precipitacion_mm)) +
  geom_point(size = 1.4) +
  scale_color_viridis_c(option = "mako", direction = -1) +
  coord_equal() +
  labs(title = paste("Observado CHIRPS -", etiqueta),
       x = "Longitud", y = "Latitud", color = "mm") +
  theme_minimal()

g_pred <- ggplot(grilla, aes(x = x, y = y, fill = precipitacion)) +
  geom_raster() +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  coord_equal() +
  labs(title = paste("Predicho kriging ordinario -", etiqueta),
       x = "Longitud", y = "Latitud", fill = "mm") +
  theme_minimal()

g_incert <- ggplot(grilla, aes(x = x, y = y, fill = desv)) +
  geom_raster() +
  scale_fill_viridis_c(option = "inferno") +
  coord_equal() +
  labs(title = "Desviacion estandar de kriging",
       x = "Longitud", y = "Latitud", fill = "mm") +
  theme_minimal()

grid.arrange(g_obs, g_pred, g_incert, ncol = 3)


# 18. Exportacion de resultados

superficie <- rast(grilla[, c("x", "y", "precipitacion")], type = "xyz",
                   crs = "EPSG:4326")
incertidumbre <- rast(grilla[, c("x", "y", "desv")], type = "xyz",
                      crs = "EPSG:4326")

writeRaster(superficie,
            paste0("prediccion_", semana_pred$anio, "_s", semana_pred$semana, ".tif"),
            overwrite = TRUE)
writeRaster(incertidumbre,
            paste0("incertidumbre_", semana_pred$anio, "_s", semana_pred$semana, ".tif"),
            overwrite = TRUE)

write.csv(resumen, "resumen_escenarios.csv", row.names = FALSE)
write.csv(resultados, "resultados_por_semana.csv", row.names = FALSE)
write.csv(d_val[, c("x", "y", "precipitacion_mm", "pred_loo", "error_loo", "var_loo")],
          "validacion_loo.csv", row.names = FALSE)

cat("\nArchivos exportados en:", getwd(), "\n")