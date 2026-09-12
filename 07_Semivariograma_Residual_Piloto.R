# =============================================================================
# 07_Semivariograma_Residual_Piloto.R
#
# Semivariograma empírico de los residuales del modelo de media del corte
# 2014-banda 38. Se comparan estructuras exponencial, esférica y gaussiana.
# El ajuste se hace sobre los residuales de sqrt(precipitación), porque esa fue
# la respuesta del mejor modelo no espacial según la validación LOO.
# =============================================================================
suppressPackageStartupMessages({ library(sp); library(gstat) })

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
                 "semivariograma_resultados")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(IN)

# Coordenadas métricas centradas.
lat0 <- mean(d$y)
d$x_km <- (d$x - mean(d$x)) * 111.32 * cos(lat0 * pi / 180)
d$y_km <- (d$y - lat0) * 111.32

# El archivo guarda el residual en escala raíz cuadrada.
spdf <- SpatialPointsDataFrame(coords = d[, c("x_km", "y_km")],
                               data = data.frame(residual_ganador = d$residual_ganador))
proj4string(spdf) <- CRS(NA_character_)

dist_max <- max(spDists(spdf))
cutoff <- 0.85 * dist_max
n_bins <- 15
vg_emp <- variogram(residual_ganador ~ 1, data = spdf,
                    cutoff = cutoff, width = cutoff / n_bins,
                    cloud = FALSE)
vg_emp <- vg_emp[is.finite(vg_emp$gamma) & vg_emp$np > 0, ]
write.csv(vg_emp[, c("dist", "gamma", "np")],
          file.path(OUT, "semivariograma_empirico_residuales.csv"),
          row.names = FALSE)

# Valores iniciales guiados por la varianza residual y la distancia máxima.
var_res <- var(d$residual_ganador)
rango_ini <- cutoff / 3
nugget_ini <- 0.10 * var_res
psill_ini <- 0.90 * var_res
modelos <- c(exponencial = "Exp", esferico = "Sph", gaussiano = "Gau")

ajustes <- lapply(modelos, function(tipo) {
  inicial <- vgm(psill = psill_ini, model = tipo, range = rango_ini,
                 nugget = nugget_ini)
  tryCatch(fit.variogram(vg_emp, inicial, fit.method = 7),
           error = function(e) NULL)
})
names(ajustes) <- names(modelos)

tabla_ajustes <- do.call(rbind, lapply(names(ajustes), function(nm) {
  a <- ajustes[[nm]]
  if (is.null(a)) return(data.frame(modelo = nm, nugget = NA, psill = NA,
    range_km = NA, SSErr = NA, stringsAsFactors = FALSE))
  data.frame(modelo = nm,
    nugget = a$psill[a$model == "Nug"][1],
    psill = a$psill[a$model != "Nug"][1],
    range_km = a$range[a$model != "Nug"][1],
    SSErr = attr(a, "SSErr"), stringsAsFactors = FALSE)
}))
tabla_ajustes <- tabla_ajustes[order(tabla_ajustes$SSErr), ]
write.csv(tabla_ajustes, file.path(OUT, "comparacion_modelos_variograma.csv"),
          row.names = FALSE)
print(tabla_ajustes, row.names = FALSE)

# Gráficos: nube conceptual (muestra) y semivariograma empírico con ajustes.
set.seed(38)
pares <- combn(nrow(d), 2)
muestra <- sample(seq_len(ncol(pares)), min(12000, ncol(pares)))
h_cloud <- sqrt((d$x_km[pares[1, muestra]] - d$x_km[pares[2, muestra]])^2 +
                (d$y_km[pares[1, muestra]] - d$y_km[pares[2, muestra]])^2)
g_cloud <- (d$residual_ganador[pares[1, muestra]] -
            d$residual_ganador[pares[2, muestra]])^2 / 2

pdf(file.path(OUT, "semivariograma_piloto_comparacion.pdf"),
    width = 11, height = 8.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))
plot(h_cloud, g_cloud, pch = 16, cex = .35, col = adjustcolor("grey30", .25),
     xlab = "Distancia entre celdas (km)", ylab = "Semivarianza",
     main = "Nube de semivariograma residual")
points(vg_emp$dist, vg_emp$gamma, pch = 19, col = "navy")
plot(vg_emp$dist, vg_emp$gamma, pch = 19, col = "navy",
     xlab = "Distancia entre celdas (km)", ylab = "Semivarianza residual",
     main = "Semivariograma empírico y modelos")
cols <- c(exponencial = "firebrick", esferico = "darkgreen", gaussiano = "darkorange")
h_line <- seq(0, max(vg_emp$dist), length.out = 200)
for (nm in names(ajustes)) if (!is.null(ajustes[[nm]])) {
  lines(h_line, variogramLine(ajustes[[nm]], dist_vector = h_line)$gamma,
        col = cols[nm], lwd = 2)
}
legend("bottomright", legend = names(cols), col = cols, lwd = 2, bty = "n")
dev.off()

cat("Semivariograma calculado con", nrow(vg_emp), "intervalos de distancia.\n")
cat("Resultados guardados en:", OUT, "\n")
cat("El siguiente paso será revisar el ajuste y decidir el modelo espacial.\n")
