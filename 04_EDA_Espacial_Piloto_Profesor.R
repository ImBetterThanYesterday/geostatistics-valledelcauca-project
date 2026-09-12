# =============================================================================
# 04_EDA_Espacial_Piloto_Profesor.R
#
# EDA espacial detallado del corte piloto: 2014 - banda 38.
# Sigue los ejemplos del profesor: distribución, mapa de posting, relaciones
# con coordenadas, tendencia preliminar, residuales y correlaciones.
# Este archivo todavía no calcula semivariogramas ni hace kriging.
# =============================================================================

suppressPackageStartupMessages(library(terra))

buscar_base <- function() {
  candidatos <- unique(normalizePath(c(getwd(), file.path(getwd(), ".."),
    file.path(getwd(), "../.."), file.path(getwd(), "../../..")), mustWork = FALSE))
  ok <- vapply(candidatos, function(p) file.exists(file.path(p,
    "avance_team_proyecto1", "reinicio_desde_cero", "cortes_espaciales_resultados",
    "tres_cortes_espaciales.rds")), logical(1))
  if (!any(ok)) stop("No se encontró tres_cortes_espaciales.rds.")
  candidatos[which(ok)[1]]
}

BASE <- buscar_base()
DIR_CORTES <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                        "cortes_espaciales_resultados")
DIR_SALIDA <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                        "eda_espacial_resultados")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

cortes <- readRDS(file.path(DIR_CORTES, "tres_cortes_espaciales.rds"))
datos <- cortes$intermedia
datos <- datos[order(datos$id_celda), ]
if (nrow(datos) != 686 || length(unique(datos$anio)) != 1 ||
    length(unique(datos$semana)) != 1) {
  stop("El corte piloto no tiene 686 celdas y un solo momento.")
}

# Coordenadas centradas en kilómetros para la tendencia.
lon0 <- mean(datos$x); lat0 <- mean(datos$y)
datos$x_km <- (datos$x - lon0) * 111.32 * cos(lat0 * pi / 180)
datos$y_km <- (datos$y - lat0) * 111.32

# Estadísticos descriptivos del corte piloto.
variables <- c("precipitacion_mm", "temperatura_bilineal_c",
               "radiacion_bilineal_mj_m2_dia", "altitud_m")
estadisticos <- do.call(rbind, lapply(variables, function(v) {
  x <- datos[[v]]
  data.frame(variable = v, n = sum(is.finite(x)), media = mean(x),
    mediana = median(x), sd = sd(x), minimo = min(x), maximo = max(x))
}))
write.csv(estadisticos, file.path(DIR_SALIDA, "estadisticos_eda_piloto.csv"),
          row.names = FALSE)

# Tendencias lineales exploratorias en las coordenadas.
modelos_tendencia <- list(
  original = lm(precipitacion_mm ~ x_km + y_km, data = datos),
  raiz = lm(sqrt(precipitacion_mm) ~ x_km + y_km, data = datos),
  log1p = lm(log1p(precipitacion_mm) ~ x_km + y_km, data = datos)
)
tabla_tendencia <- do.call(rbind, lapply(names(modelos_tendencia), function(nm) {
  m <- modelos_tendencia[[nm]]
  co <- summary(m)$coefficients
  data.frame(respuesta = nm, r2 = summary(m)$r.squared,
    r2_ajustado = summary(m)$adj.r.squared, AIC = AIC(m),
    coef_x_km = unname(co["x_km", "Estimate"]),
    coef_y_km = unname(co["y_km", "Estimate"]),
    p_x = unname(co["x_km", "Pr(>|t|)"]),
    p_y = unname(co["y_km", "Pr(>|t|)"]))
}))
write.csv(tabla_tendencia,
          file.path(DIR_SALIDA, "tendencia_espacial_preliminar.csv"),
          row.names = FALSE)

modelo_exploratorio <- modelos_tendencia$original
datos$ajuste_tendencia <- fitted(modelo_exploratorio)
datos$residual_tendencia <- residuals(modelo_exploratorio)
write.csv(datos[, c("id_celda", "x", "y", "precipitacion_mm",
                    "ajuste_tendencia", "residual_tendencia")],
          file.path(DIR_SALIDA, "residuales_tendencia_piloto.csv"),
          row.names = FALSE)

vars_cor <- c("precipitacion_mm", "altitud_m", "temperatura_bilineal_c",
              "radiacion_bilineal_mj_m2_dia", "x_km", "y_km")
write.csv(round(cor(datos[, vars_cor]), 5),
          file.path(DIR_SALIDA, "correlaciones_eda_piloto.csv"))

pdf(file.path(DIR_SALIDA, "EDA_espacial_piloto_detallado.pdf"),
    width = 11, height = 8.5)

# Página 1: distribución y mapas de posting.
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.3, 1.5))
hist(datos$precipitacion_mm, breaks = 30, col = "steelblue", border = "white",
  main = "Distribución de precipitación", xlab = "mm por banda", ylab = "Celdas")
boxplot(datos$precipitacion_mm, horizontal = TRUE, col = "steelblue",
  main = "Boxplot de precipitación", xlab = "mm por banda")
r_precip <- rast(datos[, c("x", "y", "precipitacion_mm")], type = "xyz",
                 crs = "EPSG:4326")
plot(r_precip, col = hcl.colors(40, "Blues 3", rev = TRUE),
  main = "Mapa de posting: precipitación", xlab = "Longitud", ylab = "Latitud")
r_alt <- rast(datos[, c("x", "y", "altitud_m")], type = "xyz", crs = "EPSG:4326")
plot(r_alt, col = hcl.colors(40, "Terrain 2"),
  main = "Mapa de posting: altitud", xlab = "Longitud", ylab = "Latitud")
mtext("EDA espacial — corte intermedio 2014, banda 38", outer = TRUE,
      line = -1.2, font = 2, cex = 1.2)

# Página 2: relaciones preliminares, incluyendo radiación.
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3.3, 1.5))
covs <- c("x", "y", "altitud_m", "temperatura_bilineal_c",
          "radiacion_bilineal_mj_m2_dia")
labs <- c("Longitud", "Latitud", "Altitud (m)", "Temperatura (°C)",
          "Radiación (MJ/m²/día)")
cols <- c("steelblue", "steelblue", "darkgreen", "firebrick", "darkorange")
for (j in seq_along(covs)) {
  v <- covs[j]
  plot(datos[[v]], datos$precipitacion_mm, pch = 19, cex = .6,
    col = adjustcolor(cols[j], .55), xlab = labs[j],
    ylab = "Precipitación (mm)", main = paste("Precipitación vs", labs[j]))
  abline(lm(datos$precipitacion_mm ~ datos[[v]]), lwd = 2)
  grid()
}
mtext("Relaciones preliminares: asociación, no causalidad", outer = TRUE,
      line = -1.2, font = 2, cex = 1.1)

# Página 3: tendencia lineal y residuales.
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.3, 1.5))
plot(datos$x, datos$ajuste_tendencia, pch = 19, cex = .7,
  col = hcl.colors(30, "Viridis"), xlab = "Longitud",
  ylab = "Precipitación ajustada", main = "Valores ajustados por tendencia")
r_fit <- rast(datos[, c("x", "y", "ajuste_tendencia")], type = "xyz",
              crs = "EPSG:4326")
plot(r_fit, col = hcl.colors(40, "Blues 3", rev = TRUE),
  main = "Superficie de tendencia", xlab = "Longitud", ylab = "Latitud")
plot(datos$ajuste_tendencia, datos$residual_tendencia, pch = 19, cex = .6,
  col = adjustcolor("purple", .55), xlab = "Valores ajustados",
  ylab = "Residual", main = "Residuales vs ajustados")
abline(h = 0, lty = 2, lwd = 2)
r_res <- rast(datos[, c("x", "y", "residual_tendencia")], type = "xyz",
              crs = "EPSG:4326")
plot(r_res, col = hcl.colors(40, "RdYlBu", rev = TRUE),
  main = "Mapa de residuales de tendencia",
  xlab = "Longitud", ylab = "Latitud")
mtext("La tendencia revela estructura; aún no es el modelo final", outer = TRUE,
      line = -1.2, font = 2, cex = 1.1)

# Página 4: transformaciones y correlaciones.
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.3, 1.5))
hist(sqrt(datos$precipitacion_mm), breaks = 30, col = "darkseagreen3",
  border = "white", main = "Raíz cuadrada", xlab = "sqrt(precipitación)")
hist(log1p(datos$precipitacion_mm), breaks = 30, col = "tan2",
  border = "white", main = "log1p", xlab = "log(1 + precipitación)")
mc <- cor(datos[, vars_cor])
image(seq_len(ncol(mc)), seq_len(nrow(mc)), t(mc[nrow(mc):1, ]),
  axes = FALSE, col = colorRampPalette(c("steelblue", "white", "firebrick"))(101),
  zlim = c(-1, 1), main = "Correlaciones espaciales")
labs_mc <- c("P", "Alt", "Temp", "Rad", "X", "Y")
axis(1, at = seq_len(ncol(mc)), labels = labs_mc, las = 2)
axis(2, at = seq_len(nrow(mc)), labels = rev(labs_mc), las = 2)
for (i in seq_len(nrow(mc))) for (j in seq_len(ncol(mc)))
  text(j, nrow(mc) - i + 1, sprintf("%.2f", mc[i, j]), cex = .7)
plot.new()
text(.5, .65, "La selección de respuesta y covariables\nse hará después de revisar\nla validación predictiva.", cex = 1.1)
mtext("Transformaciones y correlaciones", outer = TRUE, line = -1.2,
      font = 2, cex = 1.1)

# Página 5: climatologías como superficies espaciales de referencia.
# Para una banda fija, las climatologías cambian espacialmente entre celdas,
# pero se repiten entre años. No sustituyen la precipitación observada.
clim_vars <- c("clim_precipitacion_mm", "clim_temperatura_bilineal_c",
               "clim_radiacion_bilineal_mj_m2_dia")
clim_labs <- c("Climatología precipitación (mm)",
               "Climatología temperatura (°C)",
               "Climatología radiación (MJ/m²/día)")
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3.3, 1.5))
for (j in seq_along(clim_vars)) {
  v <- clim_vars[j]
  r <- rast(datos[, c("x", "y", v)], type = "xyz", crs = "EPSG:4326")
  pal <- if (j == 1) hcl.colors(40, "Blues 3", rev = TRUE) else
    if (j == 2) hcl.colors(40, "Inferno") else hcl.colors(40, "YlOrBr")
  plot(r, col = pal, main = clim_labs[j], xlab = "Longitud", ylab = "Latitud")
  plot(datos[[v]], datos$precipitacion_mm, pch = 19, cex = .6,
       col = adjustcolor("steelblue", .55), xlab = clim_labs[j],
       ylab = "Precipitación observada (mm)", main = "Precipitación vs climatología")
  abline(lm(datos$precipitacion_mm ~ datos[[v]]), lwd = 2)
  grid()
}
mtext("Climatologías: referencia histórica por celda y banda", outer = TRUE,
      line = -1.2, font = 2, cex = 1.1)

dev.off()
cat("EDA espacial detallado terminado para 2014-banda 38.\n")
cat("Productos guardados en:", DIR_SALIDA, "\n")
