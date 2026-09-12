# =============================================================================
# 05_Modelos_Media_Piloto.R
#
# MODELOS DE MEDIA NO ESPACIALES: CORTE 2014 - BANDA 38
#
# Se comparan modelos progresivos y tres escalas de respuesta. La validación
# principal es leave-one-out por celdas (686 observaciones). Los errores siempre
# se calculan en mm originales, para que la comparación sea interpretable.
#
# Todavía no se calcula semivariograma ni kriging.
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
DIR <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                 "cortes_espaciales_resultados")
OUT <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                 "modelos_media_resultados")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

datos <- readRDS(file.path(DIR, "tres_cortes_espaciales.rds"))$intermedia
datos <- datos[order(datos$id_celda), ]
if (nrow(datos) != 686) stop("El corte piloto debe tener 686 celdas.")

# Coordenadas en kilómetros centradas.
lon0 <- mean(datos$x); lat0 <- mean(datos$y)
datos$x_km <- (datos$x - lon0) * 111.32 * cos(lat0 * pi / 180)
datos$y_km <- (datos$y - lat0) * 111.32

# Modelos progresivos. Las climatologías se dejan como escenario exploratorio:
# clim_precipitacion usa el periodo completo y puede contener información del
# año objetivo; por eso no se usará como resultado principal sin recalcularla
# excluyendo 2014.
formulas <- list(
  M0_intercepto = precipitacion_mm ~ 1,
  M1_coordenadas = precipitacion_mm ~ x_km + y_km,
  M2_altitud = precipitacion_mm ~ altitud_m,
  M3_coord_altitud = precipitacion_mm ~ x_km + y_km + altitud_m,
  M4_ambiental = precipitacion_mm ~ x_km + y_km + altitud_m +
    temperatura_bilineal_c + radiacion_bilineal_mj_m2_dia,
  M5_climatologia = precipitacion_mm ~ x_km + y_km + altitud_m +
    clim_precipitacion_mm + clim_radiacion_bilineal_mj_m2_dia
)
respuestas <- list(
  original = function(y) y,
  raiz = function(y) sqrt(y),
  log1p = function(y) log1p(y)
)
inversas <- list(
  original = function(z) z,
  raiz = function(z) pmax(z, 0)^2,
  log1p = function(z) pmax(expm1(z), 0)
)

# Métricas en unidades originales.
metricas <- function(obs, pred) {
  e <- obs - pred
  data.frame(MAE = mean(abs(e)), RMSE = sqrt(mean(e^2)),
             sesgo = mean(pred - obs),
             R2_predictivo = 1 - sum(e^2) / sum((obs - mean(obs))^2))
}

# Validación LOO para un modelo y respuesta.
validar_loo <- function(formula, transformar, invertir, d) {
  y <- transformar(d$precipitacion_mm)
  pred <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    m <- lm(formula, data = d[-i, , drop = FALSE],
            y = y[-i])
    # Reemplaza la respuesta de la fórmula por la escala transformada.
    ff <- update(formula, y ~ .)
    m <- lm(ff, data = transform_data(d[-i, , drop = FALSE], y[-i]))
    pred[i] <- invertir(predict(m, newdata = d[i, , drop = FALSE]))
  }
  list(pred = pred, metricas = metricas(d$precipitacion_mm, pred))
}

# Pequeña función para construir una columna de respuesta transformada.
transform_data <- function(d, y) {
  d$.respuesta_transformada <- y
  d
}

resultados <- list()
predicciones <- list()
k <- 1
for (nm in names(formulas)) {
  f_base <- formulas[[nm]]
  # Sustituimos la respuesta de la fórmula por .respuesta_transformada.
  f <- update(f_base, .respuesta_transformada ~ .)
  for (resp in names(respuestas)) {
    y <- respuestas[[resp]](datos$precipitacion_mm)
    pred <- numeric(nrow(datos))
    for (i in seq_len(nrow(datos))) {
      d_train <- transform_data(datos[-i, , drop = FALSE], y[-i])
      m <- lm(f, data = d_train)
      pred[i] <- inversas[[resp]](predict(m, newdata = datos[i, , drop = FALSE]))
    }
    met <- metricas(datos$precipitacion_mm, pred)
    m_full <- lm(f, data = transform_data(datos, y))
    met$modelo <- nm
    met$respuesta <- resp
    met$n_parametros <- length(coef(m_full))
    met$R2_ajustado_entrenamiento <- summary(m_full)$adj.r.squared
    met$AIC_escala_transformada <- AIC(m_full)
    resultados[[k]] <- met
    predicciones[[paste(nm, resp, sep = "__")]] <- data.frame(
      id_celda = datos$id_celda, x = datos$x, y = datos$y,
      observado_mm = datos$precipitacion_mm, predicho_mm = pred,
      error_mm = datos$precipitacion_mm - pred,
      modelo = nm, respuesta = resp
    )
    k <- k + 1
  }
}
tabla <- do.call(rbind, resultados)
tabla <- tabla[, c("modelo", "respuesta", "n_parametros",
                   "MAE", "RMSE", "sesgo", "R2_predictivo",
                   "R2_ajustado_entrenamiento", "AIC_escala_transformada")]
tabla <- tabla[order(tabla$RMSE), ]
rownames(tabla) <- NULL
write.csv(tabla, file.path(OUT, "comparacion_modelos_media_LOO.csv"), row.names = FALSE)

# M5 se conserva para documentar cuánto mejora un predictor histórico, pero no
# se considera válido como modelo principal: clim_precipitacion_mm fue calculada
# incluyendo también el año 2014 que estamos prediciendo (fuga de información).
tabla_principal <- tabla[tabla$modelo != "M5_climatologia", , drop = FALSE]
ganador <- tabla_principal[which.min(tabla_principal$RMSE), , drop = FALSE]
write.csv(ganador, file.path(OUT, "mejor_modelo_media_LOO.csv"), row.names = FALSE)

cat("Comparación LOO (ordenada por RMSE en mm):\n")
print(tabla, row.names = FALSE)
cat("\nMejor modelo principal (sin climatología de precipitación con fuga):\n")
print(ganador, row.names = FALSE)

# Coeficientes y residuales del mejor modelo según RMSE LOO.
ganador <- ganador[1, ]
clave <- paste(ganador$modelo, ganador$respuesta, sep = "__")
pred_ganador <- predicciones[[clave]]
f_ganador <- update(formulas[[ganador$modelo]], .respuesta_transformada ~ .)
y_ganador <- respuestas[[ganador$respuesta]](datos$precipitacion_mm)
modelo_ganador <- lm(f_ganador, data = transform_data(datos, y_ganador))
datos$residual_ganador <- residuals(modelo_ganador)
datos$ajuste_ganador <- fitted(modelo_ganador)
coeficientes <- data.frame(term = names(coef(modelo_ganador)),
                           estimacion = unname(coef(modelo_ganador)),
                           p_valor = summary(modelo_ganador)$coefficients[, 4])
write.csv(coeficientes, file.path(OUT, "coeficientes_mejor_modelo.csv"),
          row.names = FALSE)
write.csv(datos[, c("id_celda", "x", "y", "precipitacion_mm",
                    "ajuste_ganador", "residual_ganador")],
          file.path(OUT, "residuales_mejor_modelo.csv"), row.names = FALSE)

# Gráficos de validación y diagnóstico del ganador.
pdf(file.path(OUT, "diagnostico_modelos_media_piloto.pdf"),
    width = 11, height = 8.5)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.3, 1.5))
plot(pred_ganador$observado_mm, pred_ganador$predicho_mm, pch = 19, cex = .6,
  col = adjustcolor("steelblue", .6), xlab = "Observado (mm)",
  ylab = "Predicho LOO (mm)", main = paste("Predicción LOO:", ganador$modelo,
                                           "|", ganador$respuesta))
abline(0, 1, lty = 2, lwd = 2)
plot(pred_ganador$predicho_mm, pred_ganador$error_mm, pch = 19, cex = .6,
  col = adjustcolor("purple", .6), xlab = "Predicho LOO (mm)",
  ylab = "Error observado - predicho (mm)", main = "Errores de validación")
abline(h = 0, lty = 2, lwd = 2)
hist(datos$residual_ganador, breaks = 30, col = "grey70", border = "white",
  main = "Residuales del ajuste completo", xlab = "Residual (mm)")
r_res <- rast(datos[, c("x", "y", "residual_ganador")], type = "xyz",
              crs = "EPSG:4326")
plot(r_res, col = hcl.colors(40, "RdYlBu", rev = TRUE),
  main = "Mapa de residuales del ganador", xlab = "Longitud", ylab = "Latitud")
dev.off()

cat("\nResultados guardados en:", OUT, "\n")
cat("El semivariograma aún no se calcula: primero revise la tabla y el PDF.\n")
