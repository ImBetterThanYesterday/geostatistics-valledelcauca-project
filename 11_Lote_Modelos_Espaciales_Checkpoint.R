# =============================================================================
# 11_Lote_Modelos_Espaciales_Checkpoint.R
#
# Algoritmo reproducible para cortes año-banda de precipitación.
# CHECKPOINT INICIAL: 2010-2011 (104 cortes). Cambiar ANIOS_PROCESAR solo
# después de revisar los resultados.
#
# Respuesta única: precipitación semanal. NO se usan climatologías.
# Los modelos candidatos evalúan coordenadas, altitud, temperatura y radiación
# por separado, antes de permitir combinaciones adicionales.
# =============================================================================

suppressPackageStartupMessages({
  library(sp)
  library(gstat)
  library(spdep)
})

ANIOS_PROCESAR <- 2010:2011
ALFA_MORAN <- 0.05
VIF_LIMITE <- 10
MEJORA_MIN_KRIGING <- 0.01
N_BINS_VARIAGRAMA <- 15
NMAX_KRIGING <- 100
MULTIPLICADOR_MAX_PREDICCION <- 3

buscar_base <- function() {
  candidatos <- unique(normalizePath(c(getwd(), file.path(getwd(), ".."),
    file.path(getwd(), "../.."), file.path(getwd(), "../../..")), mustWork = FALSE))
  ok <- vapply(candidatos, function(p) file.exists(file.path(p,
    "avance_team_proyecto1", "reinicio_desde_cero", "dataset_eda",
    "dataset_final_7_variables.rds")), logical(1))
  if (!any(ok)) stop("No se encontró dataset_final_7_variables.rds.")
  candidatos[which(ok)[1]]
}

BASE <- buscar_base()
RUTA_DATA <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                       "dataset_eda", "dataset_final_7_variables.rds")
DIR_OUT <- file.path(BASE, "avance_team_proyecto1", "reinicio_desde_cero",
                     "lote_espacial_checkpoint_2010_2011")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
datos_total <- readRDS(RUTA_DATA)

metricas <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]; e <- obs - pred
  c(n = length(obs), MAE = mean(abs(e)), RMSE = sqrt(mean(e^2)),
    sesgo = mean(pred - obs),
    R2_predictivo = 1 - sum(e^2) / sum((obs - mean(obs))^2))
}

vif_y_condicion <- function(m) {
  X <- model.matrix(m)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(X) == 0) return(list(vif_max = NA_real_, condicion = NA_real_))
  if (ncol(X) == 1) return(list(vif_max = 1, condicion = 1))
  vifs <- sapply(seq_len(ncol(X)), function(j) {
    r2 <- summary(lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared
    1 / max(1 - r2, .Machine$double.eps)
  })
  list(vif_max = max(vifs), condicion = kappa(scale(X), exact = FALSE))
}

loo_mco <- function(m, escala) {
  y <- model.response(model.frame(m))
  h <- hatvalues(m)
  e_loo <- residuals(m) / pmax(1 - h, .Machine$double.eps)
  pred_escala <- y - e_loo
  if (escala == "original") pred <- pred_escala
  if (escala == "raiz") pred <- pmax(pred_escala, 0)^2
  if (escala == "log1p") pred <- pmax(expm1(pred_escala), 0)
  pred
}

preparar_corte <- function(d) {
  lat0 <- mean(d$y)
  d$x_km <- (d$x - mean(d$x)) * 111.32 * cos(lat0 * pi / 180)
  d$y_km <- (d$y - lat0) * 111.32
  d$x2_km <- d$x_km^2
  d$y2_km <- d$y_km^2
  d$xy_km <- d$x_km * d$y_km
  d
}

modelos_rhs <- list(
  # Referencias simples
  M0_intercepto = "1",
  Mxy_coordenadas = "x_km + y_km",
  Malt_altitud = "altitud_m",
  MT_temperatura = "temperatura_bilineal_c",
  MR_radiacion = "radiacion_bilineal_mj_m2_dia",
  # Combinaciones ambientales
  MTR_temperatura_radiacion = "temperatura_bilineal_c + radiacion_bilineal_mj_m2_dia",
  MTA_temperatura_altitud = "temperatura_bilineal_c + altitud_m",
  MRA_radiacion_altitud = "radiacion_bilineal_mj_m2_dia + altitud_m",
  MTRA_ambiental_completo = "temperatura_bilineal_c + radiacion_bilineal_mj_m2_dia + altitud_m",
  # Tendencia lineal espacial más covariables
  Mxy_alt_coordenadas_altitud = "x_km + y_km + altitud_m",
  Mxy_T_coordenadas_temperatura = "x_km + y_km + temperatura_bilineal_c",
  Mxy_R_coordenadas_radiacion = "x_km + y_km + radiacion_bilineal_mj_m2_dia",
  Mxy_TR_coordenadas_temp_rad = "x_km + y_km + temperatura_bilineal_c + radiacion_bilineal_mj_m2_dia",
  Mxy_TA_coordenadas_temp_alt = "x_km + y_km + temperatura_bilineal_c + altitud_m",
  Mxy_RA_coordenadas_rad_alt = "x_km + y_km + radiacion_bilineal_mj_m2_dia + altitud_m",
  Mxy_TRA_coordenadas_ambiental = "x_km + y_km + temperatura_bilineal_c + radiacion_bilineal_mj_m2_dia + altitud_m"
)

ajustar_variogramas <- function(residuales, d) {
  spdf <- SpatialPointsDataFrame(d[, c("x_km", "y_km")],
                                  data.frame(residual = residuales))
  dist_max <- max(spDists(spdf))
  cutoff <- .85 * dist_max
  vg <- variogram(residual ~ 1, spdf, cutoff = cutoff,
                  width = cutoff / N_BINS_VARIAGRAMA)
  vg <- vg[is.finite(vg$gamma) & vg$np > 0, ]
  if (nrow(vg) < 8 || max(vg$gamma) <= 0) {
    return(list(ok = FALSE, razon = "semivariograma_insuficiente", vg = vg))
  }
  vr <- var(residuales)
  tipos <- c(exponencial = "Exp", esferico = "Sph", gaussiano = "Gau")
  ajustes <- lapply(tipos, function(tp) suppressWarnings(tryCatch(
    fit.variogram(vg, vgm(.9 * vr, tp, cutoff / 3, .1 * vr), fit.method = 7),
    error = function(e) NULL)))
  names(ajustes) <- names(tipos)

  tabla <- do.call(rbind, lapply(names(ajustes), function(nm) {
    a <- ajustes[[nm]]
    if (is.null(a) || !is.finite(attr(a, "SSErr"))) {
      return(data.frame(modelo_variograma = nm, nugget = NA_real_,
        psill = NA_real_, range_km = NA_real_, SSErr = NA_real_))
    }
    data.frame(modelo_variograma = nm,
      nugget = a$psill[a$model == "Nug"][1],
      psill = a$psill[a$model != "Nug"][1],
      range_km = a$range[a$model != "Nug"][1],
      SSErr = attr(a, "SSErr"))
  }))
  tabla_ok <- tabla[is.finite(tabla$SSErr) & tabla$nugget >= 0 &
                      tabla$psill > 0 & tabla$range_km > 0, , drop = FALSE]

  if (nrow(tabla_ok) == 0) return(list(ok = FALSE, razon = "sin_convergencia", vg = vg, tabla = tabla))

  mejor <- tabla_ok[which.min(tabla_ok$SSErr), , drop = FALSE]
  list(ok = TRUE, vg = vg, tabla = tabla, mejor = mejor,
       ajuste = ajustes[[mejor$modelo_variograma]])
}

diagnosticar_anisotropia <- function(residuales, d) {
  spdf <- SpatialPointsDataFrame(d[, c("x_km", "y_km")],
                                  data.frame(residual = residuales))
  vm <- variogram(residual ~ 1, spdf, alpha = c(0,45,90,135),
                  tol.hor = 22.5, cutoff = .6 * max(spDists(spdf)), width = 10)
  # Solo se comparan retardos cortos comunes. Comparar la mediana de todos
  # los retardos mezcla distancias distintas según dirección y sobre-detecta
  # anisotropía en dominios de forma irregular como Valle del Cauca.
  limite_corto <- median(vm$dist, na.rm = TRUE)
  med <- aggregate(gamma ~ dir.hor, vm[vm$dist <= limite_corto, ], median)
  ratio <- if (nrow(med) >= 2 && min(med$gamma) > 0) max(med$gamma) / min(med$gamma) else NA_real_
  # Es una alerta para inspección posterior; en este checkpoint no se ajusta
  # anisotropía automáticamente. Primero se compara el modelo isotrópico.
  list(ratio = ratio, bandera = is.finite(ratio) && ratio > 2)
}

procesar_corte <- function(d, anio, banda) {
  salida <- data.frame(anio = anio, banda = banda, n_celdas = nrow(d),
                       estado = "INICIADO", motivo = "", stringsAsFactors = FALSE)
  requeridas <- c("id_celda","x","y","altitud_m","precipitacion_mm",
                  "temperatura_bilineal_c","radiacion_bilineal_mj_m2_dia")
  if (nrow(d) != 686 || length(unique(d$id_celda)) != 686 ||
      anyDuplicated(d[, c("x","y")]) || anyNA(d[, requeridas]) ||
      any(d$precipitacion_mm < 0)) {
    salida$estado <- "FALLIDO_CONTROL"; salida$motivo <- "estructura_o_valores_invalidos"
    return(list(resumen = salida))
  }
  d <- preparar_corte(d)
  escalas <- list(original = d$precipitacion_mm,
                  raiz = sqrt(d$precipitacion_mm),
                  log1p = log1p(d$precipitacion_mm))
  candidatos <- list(); cc <- 1
  for (escala in names(escalas)) {
    d$.respuesta <- escalas[[escala]]
    for (nm in names(modelos_rhs)) {
      f <- as.formula(paste(".respuesta ~", modelos_rhs[[nm]]))
      m <- tryCatch(lm(f, data = d), error = function(e) NULL)
      if (is.null(m)) next
      vc <- vif_y_condicion(m)
      pred <- loo_mco(m, escala)
      mt <- metricas(d$precipitacion_mm, pred)
      candidatos[[cc]] <- data.frame(anio = anio, banda = banda, escala = escala,
        modelo_media = nm, formula = deparse(f), vif_max = vc$vif_max,
        condicion = vc$condicion, colineal = is.finite(vc$vif_max) && vc$vif_max > VIF_LIMITE,
        t(mt), stringsAsFactors = FALSE)
      cc <- cc + 1
    }
  }
  tabla_candidatos <- do.call(rbind, candidatos)
  validos <- tabla_candidatos[!tabla_candidatos$colineal, , drop = FALSE]
  if (nrow(validos) == 0) {
    salida$estado <- "FALLIDO_COLINEALIDAD"; salida$motivo <- "todos_los_modelos_vif_alto"
    return(list(resumen = salida, candidatos = tabla_candidatos))
  }
  # M0 original es la referencia común de mejora.
  base <- tabla_candidatos[tabla_candidatos$modelo_media == "M0_intercepto" &
                            tabla_candidatos$escala == "original", ]
  ganador <- validos[which.min(validos$RMSE), , drop = FALSE]
  salida <- cbind(salida, ganador[, c("escala","modelo_media","formula","vif_max",
    "condicion","MAE","RMSE","sesgo","R2_predictivo")])
  salida$mejora_mco_vs_base <- if (nrow(base)) 1 - salida$RMSE / base$RMSE else NA_real_

  d$.respuesta <- escalas[[ganador$escala]]
  f_ganador <- as.formula(ganador$formula)
  m_ganador <- lm(f_ganador, data = d)
  residuos <- residuals(m_ganador)
  moran_res <- tryCatch({
    nb <- knn2nb(knearneigh(as.matrix(d[,c("x_km","y_km")]), k = 4))
    mt <- moran.test(residuos, nb2listw(nb, style="W"), alternative="greater")
    c(I = unname(mt$estimate[["Moran I statistic"]]), p = mt$p.value)
  }, error = function(e) c(I = NA_real_, p = NA_real_))
  salida$moran_I <- moran_res["I"]; salida$moran_p <- moran_res["p"]
  salida$residual_sd <- sd(residuos)

  if (!is.finite(salida$moran_p) || salida$moran_p >= ALFA_MORAN) {
    salida$estado <- "FINAL_MCO"; salida$motivo <- "sin_autocorrelacion_residual"
    return(list(resumen = salida, candidatos = tabla_candidatos))
  }

  aniso <- diagnosticar_anisotropia(residuos, d)
  salida$anisotropia_ratio <- aniso$ratio
  salida$anisotropia_bandera <- aniso$bandera
  vg_out <- ajustar_variogramas(residuos, d)
  if (!vg_out$ok) {
    salida$estado <- "REVISAR_VARIAGRAMA"; salida$motivo <- vg_out$razon
    return(list(resumen = salida, candidatos = tabla_candidatos,
                variograma = vg_out$tabla))
  }
  salida$modelo_variograma <- vg_out$mejor$modelo_variograma
  salida$nugget <- vg_out$mejor$nugget
  salida$psill <- vg_out$mejor$psill
  salida$range_km <- vg_out$mejor$range_km
  salida$SSErr_variograma <- vg_out$mejor$SSErr

  # CV de kriging con variograma fijo en el corte; se etiqueta explícitamente.
  columnas_sp <- c(".respuesta","x_km","y_km","x2_km","y2_km","xy_km",
                   "altitud_m","temperatura_bilineal_c","radiacion_bilineal_mj_m2_dia")
  spdf <- SpatialPointsDataFrame(d[,c("x_km","y_km")], data=d[,columnas_sp])
  tipo_kriging <- if (ganador$modelo_media == "M0_intercepto") "ordinario" else "universal"
  f_kriging <- if (tipo_kriging == "ordinario") .respuesta ~ 1 else f_ganador
  cv <- tryCatch(krige.cv(f_kriging, spdf, model = vg_out$ajuste,
                           nfold = nrow(d), nmax = NMAX_KRIGING, debug.level = 0),
                 error = function(e) NULL)
  if (is.null(cv)) {
    salida$estado <- "REVISAR_KRIGING"; salida$motivo <- "error_krige_cv"
    return(list(resumen = salida, candidatos = tabla_candidatos,
                variograma = vg_out$tabla))
  }
  pred_esc <- cv$var1.pred
  prop_pred_neg_escala <- mean(pred_esc < 0, na.rm = TRUE)
  pred <- if (ganador$escala == "original") pred_esc else
    if (ganador$escala == "raiz") pmax(pred_esc,0)^2 else pmax(expm1(pred_esc),0)
  mk <- metricas(d$precipitacion_mm, pred)
  salida$tipo_kriging <- tipo_kriging
  salida$MAE_kriging <- mk["MAE"]; salida$RMSE_kriging <- mk["RMSE"]
  salida$sesgo_kriging <- mk["sesgo"]; salida$R2_kriging <- mk["R2_predictivo"]
  salida$n_predicciones_kriging_validas <- mk["n"]
  salida$max_prediccion_kriging_mm <- if (any(is.finite(pred))) max(pred, na.rm = TRUE) else NA_real_
  salida$proporcion_pred_negativa_escala <- prop_pred_neg_escala
  salida$cobertura_95_escala <- mean(abs(cv$residual) <= 1.96 * sqrt(pmax(cv$var1.var,0)))
  salida$validacion_kriging <- "LOO_con_variograma_fijo"
  # Predicciones no finitas o extremadamente mayores que los datos del corte
  # son inestabilidad numérica, no evidencia a favor de kriging. Se conserva
  # el MCO y se registra el caso para inspección.
  limite_pred <- MULTIPLICADOR_MAX_PREDICCION * max(d$precipitacion_mm)
  if (salida$n_predicciones_kriging_validas < nrow(d) ||
      !is.finite(salida$max_prediccion_kriging_mm) ||
      salida$max_prediccion_kriging_mm > limite_pred) {
    salida$estado <- "FINAL_MCO_KRIGING_INESTABLE"
    salida$motivo <- "predicciones_kriging_no_finitas_o_extremas"
  } else if (is.finite(salida$RMSE_kriging) &&
      salida$RMSE_kriging < salida$RMSE * (1 - MEJORA_MIN_KRIGING)) {
    salida$estado <- if (isTRUE(aniso$bandera)) "FINAL_KRIGING_REVISAR_ANISOTROPIA" else "FINAL_KRIGING"
    salida$motivo <- "kriging_mejora_rmse_mco_mas_1pct"
  } else {
    salida$estado <- "FINAL_MCO"
    salida$motivo <- "kriging_no_mejora_rmse_mco_mas_1pct"
  }
  list(resumen = salida, candidatos = tabla_candidatos, variograma = vg_out$tabla)
}

cortes <- expand.grid(anio = ANIOS_PROCESAR, banda = 1:52)
resumenes <- list(); candidatos_todos <- list(); variogramas_todos <- list()
for (i in seq_len(nrow(cortes))) {
  a <- cortes$anio[i]; b <- cortes$banda[i]
  cat(sprintf("[%03d/%03d] %d banda %02d\n", i, nrow(cortes), a, b))
  d <- datos_total[datos_total$anio == a & datos_total$semana == b, ]
  ans <- tryCatch(procesar_corte(d, a, b), error = function(e)
    list(resumen = data.frame(anio=a,banda=b,n_celdas=nrow(d),estado="ERROR",
      motivo=conditionMessage(e),stringsAsFactors=FALSE)))
  resumenes[[i]] <- ans$resumen
  if (!is.null(ans$candidatos)) candidatos_todos[[i]] <- ans$candidatos
  if (!is.null(ans$variograma)) {
    ans$variograma$anio <- a; ans$variograma$banda <- b
    variogramas_todos[[i]] <- ans$variograma
  }
  if (i %% 10 == 0) {
    write.csv(do.call(rbind,resumenes), file.path(DIR_OUT,"checkpoint_resumen_parcial.csv"),
              row.names=FALSE)
  }
}

resumen_final <- do.call(rbind,resumenes)
write.csv(resumen_final,file.path(DIR_OUT,"resumen_decisiones_2010_2011.csv"),row.names=FALSE)
if(length(candidatos_todos)) write.csv(do.call(rbind,candidatos_todos),
  file.path(DIR_OUT,"modelos_media_candidatos_2010_2011.csv"),row.names=FALSE)
if(length(variogramas_todos)) write.csv(do.call(rbind,variogramas_todos),
  file.path(DIR_OUT,"variogramas_candidatos_2010_2011.csv"),row.names=FALSE)

sink(file.path(DIR_OUT,"CONCLUSION_CHECKPOINT.txt"))
cat("CHECKPOINT: 2010-2011\n\n")
print(table(resumen_final$estado,useNA="ifany"))
cat("\nEl lote completo 2010-2025 NO se ejecuta hasta revisar este resumen.\n")
sink()

cat("\nCheckpoint terminado. Revisar:", DIR_OUT, "\n")
print(table(resumen_final$estado,useNA="ifany"))
