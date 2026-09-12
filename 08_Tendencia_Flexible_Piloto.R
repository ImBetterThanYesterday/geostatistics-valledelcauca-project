# =============================================================================
# 08_Tendencia_Flexible_Piloto.R
#
# Compara la tendencia espacial lineal con una tendencia cuadrática en el corte
# 2014-banda 38. No hace kriging todavía.
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(spdep); library(gstat); library(sp) })

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
D <- readRDS(file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                       "cortes_espaciales_resultados", "tres_cortes_espaciales.rds"))$intermedia
D <- D[order(D$id_celda), ]
lat0 <- mean(D$y)
D$x_km <- (D$x - mean(D$x)) * 111.32 * cos(lat0*pi/180)
D$y_km <- (D$y - lat0) * 111.32
D$x2_km <- D$x_km^2
D$y2_km <- D$y_km^2
D$xy_km <- D$x_km * D$y_km
D$.y <- sqrt(D$precipitacion_mm)

OUT <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                 "tendencia_flexible_resultados")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

forms <- list(
  lineal_ambiental = .y ~ x_km + y_km + altitud_m +
    temperatura_bilineal_c + radiacion_bilineal_mj_m2_dia,
  cuadratica_ambiental = .y ~ x_km + y_km + x2_km + y2_km + xy_km +
    altitud_m + temperatura_bilineal_c + radiacion_bilineal_mj_m2_dia
)

metricas <- function(obs, pred) {
  e <- obs - pred
  c(MAE = mean(abs(e)), RMSE = sqrt(mean(e^2)),
    sesgo = mean(pred - obs),
    R2_predictivo = 1 - sum(e^2)/sum((obs-mean(obs))^2))
}
preds <- list(); rows <- list()
for (nm in names(forms)) {
  f <- forms[[nm]]
  pred <- numeric(nrow(D))
  for (i in seq_len(nrow(D))) {
    m <- lm(f, data = D[-i, , drop = FALSE])
    pred[i] <- pmax(predict(m, newdata = D[i, , drop = FALSE]), 0)^2
  }
  z <- metricas(D$precipitacion_mm, pred)
  rows[[nm]] <- data.frame(modelo = nm, t(z), r2_ajustado = summary(lm(f,D))$adj.r.squared)
  preds[[nm]] <- pred
}
comparacion <- do.call(rbind, rows)
write.csv(comparacion, file.path(OUT, "comparacion_tendencia_lineal_cuadratica_LOO.csv"),
          row.names = FALSE)
print(comparacion, row.names = FALSE)

# Residuales transformados y Moran para cada tendencia.
moran_rows <- list(); residuales <- list()
for (nm in names(forms)) {
  m <- lm(forms[[nm]], data = D)
  r <- residuals(m); residuales[[nm]] <- r
  xy <- as.matrix(D[, c("x_km","y_km")])
  mor <- lapply(c(4,8), function(k) {
    nb <- knn2nb(knearneigh(xy, k=k))
    lw <- nb2listw(nb, style="W", zero.policy=TRUE)
    mt <- moran.test(r, lw, zero.policy=TRUE, alternative="greater")
    data.frame(modelo=nm, vecinos=k, moran_I=unname(mt$estimate[["Moran I statistic"]]),
               p_valor=mt$p.value)
  })
  moran_rows <- c(moran_rows, mor)
  write.csv(data.frame(id_celda=D$id_celda,x=D$x,y=D$y,residual=r),
            file.path(OUT, paste0("residuales_",nm,".csv")), row.names=FALSE)
}
moran_tabla <- do.call(rbind, moran_rows)
write.csv(moran_tabla, file.path(OUT, "moran_tendencias.csv"), row.names=FALSE)
print(moran_tabla, row.names=FALSE)

# Semivariograma empírico de cada tendencia para comparar su forma.
vg_rows <- list()
for (nm in names(residuales)) {
  spdf <- SpatialPointsDataFrame(coords=D[,c("x_km","y_km")],
                                  data=data.frame(residual=residuales[[nm]]))
  vg <- variogram(residual ~ 1, spdf, cutoff=0.85*max(spDists(spdf)),
                  width=(0.85*max(spDists(spdf)))/15)
  vg$modelo <- nm
  vg_rows[[nm]] <- vg[,c("dist","gamma","np","modelo")]
}
vg_tabla <- do.call(rbind, vg_rows)
write.csv(vg_tabla, file.path(OUT, "semivariogramas_comparados.csv"), row.names=FALSE)

pdf(file.path(OUT, "comparacion_tendencias_residuales.pdf"), width=11, height=8.5)
par(mfrow=c(2,2), mar=c(4.5,4.5,3.3,1.5))
for (nm in names(residuales)) {
  r <- residuales[[nm]]
  plot(D$ajuste_dummy <- fitted(lm(forms[[nm]],D)), r, pch=19, cex=.55,
       col=adjustcolor("purple",.5), xlab="Ajuste en escala raíz",
       ylab="Residual", main=paste("Residuales:",nm))
  abline(h=0,lty=2,lwd=2)
  spdf <- SpatialPointsDataFrame(coords=D[,c("x_km","y_km")],
                                  data=data.frame(residual=r))
  plot(rast(data.frame(x=D$x,y=D$y,residual=r),type="xyz",crs="EPSG:4326"),
       col=hcl.colors(40,"RdYlBu",rev=TRUE), main=paste("Mapa:",nm),
       xlab="Longitud",ylab="Latitud")
}
dev.off()
cat("Resultados guardados en:", OUT, "\n")

