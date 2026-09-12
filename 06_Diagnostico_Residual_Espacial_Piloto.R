# =============================================================================
# 06_Diagnostico_Residual_Espacial_Piloto.R
#
# Prueba de autocorrelación espacial de los residuales del mejor modelo de
# media del corte 2014-banda 38. No ajusta todavía un variograma ni kriging.
# Se evalúan vecindades de 4 y 8 celdas cercanas como análisis de sensibilidad.
# =============================================================================
suppressPackageStartupMessages(library(spdep))

buscar_base <- function() {
  candidatos <- unique(normalizePath(c(getwd(), file.path(getwd(), ".."),
    file.path(getwd(), "../.."), file.path(getwd(), "../../..")), mustWork = FALSE))
  ok <- vapply(candidatos, function(p) file.exists(file.path(p,
    "avance_team_proyecto1", "reinicio_desde_cero", "modelos_media_resultados",
    "residuales_mejor_modelo.csv")), logical(1))
  if (!any(ok)) stop("No se encontró residuales_mejor_modelo.csv.")
  candidatos[which(ok)[1]]
}
BASE <- buscar_base()
IN <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                "modelos_media_resultados", "residuales_mejor_modelo.csv")
OUT <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                 "dependencia_residual_resultados")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(IN)

# Distancias en km para que el umbral sea interpretable.
lat0 <- mean(d$y)
d$x_km <- (d$x - mean(d$x)) * 111.32 * cos(lat0 * pi / 180)
d$y_km <- (d$y - lat0) * 111.32

# k=4 aproxima vecinos ortogonales; k=8 incluye diagonales.
resultados <- do.call(rbind, lapply(c(4, 8), function(k) {
  kn <- knn2nb(knearneigh(as.matrix(d[, c("x_km", "y_km")]), k = k))
  lw <- nb2listw(kn, style = "W", zero.policy = TRUE)
  mr <- moran.test(d$residual_ganador, lw, zero.policy = TRUE,
                   alternative = "greater", randomisation = TRUE)
  mp <- moran.test(d$precipitacion_mm, lw, zero.policy = TRUE,
                   alternative = "greater", randomisation = TRUE)
  data.frame(vecinos = k, variable = c("residual_ganador", "precipitacion_mm"),
    moran_I = c(unname(mr$estimate[["Moran I statistic"]]),
                unname(mp$estimate[["Moran I statistic"]])),
    esperanza_I = c(unname(mr$estimate[["Expectation"]]),
                    unname(mp$estimate[["Expectation"]])),
    p_valor = c(mr$p.value, mp$p.value),
    n_celdas = nrow(d))
}))
write.csv(resultados, file.path(OUT, "moran_residual_vs_crudo.csv"), row.names = FALSE)
print(resultados, row.names = FALSE)
cat("\nConclusión: si Moran de residuales sigue siendo positivo y significativo,\n",
    "queda dependencia espacial y corresponde construir el semivariograma.\n", sep = "")

