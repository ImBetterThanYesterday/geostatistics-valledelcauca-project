# =============================================================================
# 13_Ver_Graficas_Resultados_Lote.R
# Visualizacion del lote espacial (no reentrena ningun modelo).
#
# Uso desde RStudio, con el proyecto abierto:
# source("avance_team_proyecto1/reinicio_desde_cero/13_Ver_Graficas_Resultados_Lote.R")
#
# Cambie CORTE para inspeccionar cualquier anio y semana que haya terminado
# con estado FINAL_KRIGING.
# =============================================================================

buscar_raiz <- function() {
  candidatos <- unique(normalizePath(c(getwd(), file.path(getwd(), ".."),
                                       file.path(getwd(), "../..")), mustWork = FALSE))
  ok <- vapply(candidatos, function(p) {
    file.exists(file.path(p, "avance_team_proyecto1", "reinicio_desde_cero",
                          "resultados_baseR_c2", "resumen_decisiones.csv"))
  }, logical(1))
  if (!any(ok)) stop("No se encontro la raiz del proyecto. Abra la carpeta Proyecto1 en RStudio.")
  candidatos[which(ok)[1]]
}

RAIZ_PROYECTO <- buscar_raiz()
RUTA_RESULTADOS <- file.path(RAIZ_PROYECTO, "avance_team_proyecto1", "reinicio_desde_cero",
                              "resultados_baseR_c2")
RUTA_GRAFICOS <- file.path(RUTA_RESULTADOS, "graficos_resumen")
dir.create(RUTA_GRAFICOS, recursive = TRUE, showWarnings = FALSE)

# Corte de ejemplo. Puede cambiarlo, por ejemplo: CORTE <- c(2022, 10)
CORTE <- c(2014L, 38L)

resumen <- read.csv(file.path(RUTA_RESULTADOS, "resumen_decisiones.csv"), check.names = FALSE)
predicciones <- readRDS(file.path(RUTA_RESULTADOS, "predicciones.rds"))

guardar_png <- function(nombre, expr, ancho = 1800, alto = 1200) {
  png(file.path(RUTA_GRAFICOS, nombre), width = ancho, height = alto, res = 160)
  on.exit(dev.off(), add = TRUE)
  eval.parent(substitute(expr))
}

paleta_mapa <- function(x, n = 80, esquema = c("azul", "residual")) {
  esquema <- match.arg(esquema)
  if (esquema == "azul") {
    br <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1)
    list(ind = pmax(1, pmin(n, findInterval(x, br, all.inside = TRUE))),
         col = hcl.colors(n, "Blues 3", rev = TRUE), br = br)
  } else {
    lim <- max(abs(x), na.rm = TRUE)
    br <- seq(-lim, lim, length.out = n + 1)
    list(ind = pmax(1, pmin(n, findInterval(x, br, all.inside = TRUE))),
         col = hcl.colors(n, "Blue-Red 3"), br = br)
  }
}

finales <- resumen[resumen$estado == "FINAL_KRIGING", ]
cat("\nResultados disponibles:", nrow(resumen), "cortes;", nrow(finales), "con kriging final.\n")
cat("Las figuras se guardaran en:\n", RUTA_GRAFICOS, "\n\n", sep = "")

# 1. Estado del lote
guardar_png("01_estado_lote.png", {
  conteo <- table(resumen$estado)
  barplot(conteo, col = c("steelblue", "darkorange")[seq_along(conteo)],
          main = "Estado de los 780 cortes anio × semana", ylab = "Numero de cortes")
  text(seq_along(conteo), conteo, labels = conteo, pos = 3)
})

# 2. Que modelos y escalas fueron seleccionados
guardar_png("02_modelos_y_transformaciones.png", {
  par(mfrow = c(1, 2), mar = c(10, 5, 4, 1))
  tm <- sort(table(finales$modelo), decreasing = TRUE)
  barplot(tm, las = 2, col = "steelblue", main = "Modelo de media elegido", ylab = "Cortes")
  te <- table(finales$escala)
  barplot(te, col = c("grey55", "seagreen3", "mediumpurple3"),
          main = "Escala de respuesta elegida", ylab = "Cortes")
  text(seq_along(te), te, labels = te, pos = 3)
})

# 3. Semivariograma elegido y sus parametros
guardar_png("03_semivariogramas.png", {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
  tf <- table(finales$familia)
  barplot(tf, col = c("tomato3", "goldenrod2", "mediumpurple3"),
          main = "Familia de semivariograma elegida", ylab = "Cortes")
  boxplot(range_km ~ familia, data = finales, col = c("tomato3", "goldenrod2", "mediumpurple3"),
          main = "Alcance espacial por familia", xlab = "Familia", ylab = "Alcance (km)")
})

# 4. Metricas predictivas: MCO vs MCG/kriging
guardar_png("04_metricas_mco_vs_kriging.png", {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
  ok <- is.finite(finales$RMSE) & is.finite(finales$RMSE_mcg_kriging)
  boxplot(finales$RMSE[ok], finales$RMSE_mcg_kriging[ok], names = c("MCO", "MCG + kriging"),
          col = c("grey70", "steelblue"), main = "RMSE por corte", ylab = "RMSE (mm)")
  mejora <- finales$RMSE[ok] - finales$RMSE_mcg_kriging[ok]
  hist(mejora, breaks = 35, col = "steelblue", border = "white",
       main = "Mejora de RMSE al usar kriging", xlab = "RMSE MCO − RMSE kriging (mm)")
  abline(v = 0, lty = 2, lwd = 2)
})

# 5. Moran de residuales del modelo de media
guardar_png("05_dependencia_espacial_residual.png", {
  hist(finales$moran_I, breaks = 35, col = "mediumpurple3", border = "white",
       main = "Moran I de residuales MCO", xlab = "Moran I")
  mtext("Valores positivos indican agrupamiento espacial residual", side = 3, line = 0.3, cex = 0.85)
})

# 6. Un corte completo: observado, prediccion, residual y comparacion
sub <- predicciones[predicciones$anio == CORTE[1] & predicciones$semana == CORTE[2], ]
if (nrow(sub) == 0) {
  warning("El corte solicitado no tiene predicciones. Elija uno con FINAL_KRIGING en resumen_decisiones.csv.")
} else {
  sub$residual_mm <- sub$observado_mm - sub$prediccion_mm
  po <- paleta_mapa(sub$observado_mm, esquema = "azul")
  pp <- paleta_mapa(sub$prediccion_mm, esquema = "azul")
  pr <- paleta_mapa(sub$residual_mm, esquema = "residual")
  guardar_png(sprintf("06_corte_%d_semana_%02d.png", CORTE[1], CORTE[2]), {
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    plot(sub$x, sub$y, pch = 15, cex = 1.25, col = po$col[po$ind], asp = 1,
         xlab = "Longitud", ylab = "Latitud", main = "Precipitacion observada")
    plot(sub$x, sub$y, pch = 15, cex = 1.25, col = pp$col[pp$ind], asp = 1,
         xlab = "Longitud", ylab = "Latitud", main = "Prediccion por kriging")
    plot(sub$x, sub$y, pch = 15, cex = 1.25, col = pr$col[pr$ind], asp = 1,
         xlab = "Longitud", ylab = "Latitud", main = "Residual: observado − predicho")
    plot(sub$observado_mm, sub$prediccion_mm, pch = 16, col = rgb(0.1, 0.3, 0.7, 0.5),
         xlab = "Precipitacion observada (mm)", ylab = "Prediccion LOO (mm)",
         main = "Validacion punto a punto")
    abline(0, 1, lty = 2, lwd = 2, col = "red3")
    mtext(sprintf("Corte espacial: %d, semana ISO %02d", CORTE[1], CORTE[2]), outer = TRUE, line = -1.2)
  }, ancho = 1800, alto = 1500)
}

cat("Listo. Abra la pestana Files y entre a resultados_baseR_c2/graficos_resumen.\n")
