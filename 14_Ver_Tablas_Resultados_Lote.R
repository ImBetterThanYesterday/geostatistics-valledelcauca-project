# =============================================================================
# 14_Ver_Tablas_Resultados_Lote.R
# Tablas legibles del lote espacial. No reentrena ningun modelo.
# =============================================================================

buscar_raiz <- function() {
  candidatos <- unique(normalizePath(c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../..")), mustWork = FALSE))
  ok <- vapply(candidatos, function(p) file.exists(file.path(
    p, "avance_team_proyecto1", "reinicio_desde_cero", "resultados_baseR_c2", "resumen_decisiones.csv")), logical(1))
  if (!any(ok)) stop("No se encontro la raiz del proyecto. Abra Proyecto1 en RStudio.")
  candidatos[which(ok)[1]]
}

# Cambie este corte para ver su competencia completa entre los 48 modelos.
CORTE <- c(2014L, 38L)

RAIZ_PROYECTO <- buscar_raiz()
RUTA <- file.path(RAIZ_PROYECTO, "avance_team_proyecto1", "reinicio_desde_cero", "resultados_baseR_c2")
SALIDA <- file.path(RUTA, "tablas_resumen")
dir.create(SALIDA, recursive = TRUE, showWarnings = FALSE)

resumen <- read.csv(file.path(RUTA, "resumen_decisiones.csv"), check.names = FALSE)
candidatos <- read.csv(file.path(RUTA, "modelos_candidatos.csv"), check.names = FALSE)
finales <- resumen[resumen$estado == "FINAL_KRIGING", ]

# Tabla A: una fila por corte. Esta es la tabla principal: el ganador ya esta
# seleccionado por menor RMSE LOO del modelo de media y luego se compara con
# kriging bajo la misma validacion.
tabla_cortes <- resumen[order(resumen$anio, resumen$semana), c(
  "anio", "semana", "estado", "modelo", "escala", "formula", "vif_max",
  "MAE", "RMSE", "R2", "moran_I", "moran_p", "familia", "nugget", "psill",
  "range_km", "tipo_kriging", "MAE_mcg_kriging", "RMSE_mcg_kriging",
  "R2_predictivo_mcg_kriging", "cobertura95_escala"
)]
names(tabla_cortes)[names(tabla_cortes) == "modelo"] <- "modelo_medio_ganador"
names(tabla_cortes)[names(tabla_cortes) == "escala"] <- "transformacion_ganadora"
write.csv(tabla_cortes, file.path(SALIDA, "01_ganador_por_corte.csv"), row.names = FALSE)

# Tabla B: frecuencia de cada modelo ganador y rendimiento promedio. No sirve
# para reemplazar la tabla por corte: resume los 769 resultados aceptados.
agregar_modelo <- function(z) {
  data.frame(
    modelo_medio_ganador = z$modelo[1], cortes = nrow(z),
    porcentaje = 100 * nrow(z) / nrow(finales),
    MAE_medio_MCO = mean(z$MAE, na.rm = TRUE),
    RMSE_medio_MCO = mean(z$RMSE, na.rm = TRUE),
    MAE_medio_kriging = mean(z$MAE_mcg_kriging, na.rm = TRUE),
    RMSE_medio_kriging = mean(z$RMSE_mcg_kriging, na.rm = TRUE),
    mejora_RMSE_media_mm = mean(z$RMSE - z$RMSE_mcg_kriging, na.rm = TRUE)
  )
}
tabla_modelos <- do.call(rbind, lapply(split(finales, finales$modelo), agregar_modelo))
tabla_modelos <- tabla_modelos[order(-tabla_modelos$cortes), ]
write.csv(tabla_modelos, file.path(SALIDA, "02_resumen_modelos_ganadores.csv"), row.names = FALSE)

tabla_escala <- as.data.frame(table(finales$escala), stringsAsFactors = FALSE)
names(tabla_escala) <- c("transformacion_ganadora", "cortes")
tabla_escala$porcentaje <- 100 * tabla_escala$cortes / nrow(finales)
write.csv(tabla_escala, file.path(SALIDA, "03_resumen_transformaciones.csv"), row.names = FALSE)

tabla_variograma <- as.data.frame(table(finales$familia), stringsAsFactors = FALSE)
names(tabla_variograma) <- c("semivariograma_ganador", "cortes")
tabla_variograma$porcentaje <- 100 * tabla_variograma$cortes / nrow(finales)
write.csv(tabla_variograma, file.path(SALIDA, "04_resumen_semivariogramas.csv"), row.names = FALSE)

# Tabla E: el duelo de candidatos para UN corte. Menor RMSE es el ganador;
# colineal=TRUE indica que fue excluido antes de competir.
tabla_competencia <- candidatos[candidatos$anio == CORTE[1] & candidatos$semana == CORTE[2], ]
tabla_competencia$estado_candidato <- ifelse(tabla_competencia$colineal, "EXCLUIDO: VIF > 10", "VALIDO")
tabla_competencia <- tabla_competencia[order(tabla_competencia$colineal, tabla_competencia$RMSE), ]
write.csv(tabla_competencia, file.path(SALIDA, sprintf("05_competencia_%d_semana_%02d.csv", CORTE[1], CORTE[2])), row.names = FALSE)

cat("\n================ GANADORES: VISTA GLOBAL ================\n")
cat("Un corte = un anio y una semana ISO.\n")
cat("El modelo de media ganador minimiza RMSE LOO entre candidatos sin VIF > 10.\n")
cat("El metodo final solo es kriging si mejora ese RMSE en el mismo corte.\n\n")
cat("Estados finales:\n"); print(table(resumen$estado))
cat("\nModelos de media ganadores:\n"); print(tabla_modelos, row.names = FALSE)
cat("\nTransformaciones ganadoras:\n"); print(tabla_escala, row.names = FALSE)
cat("\nSemivariogramas ganadores:\n"); print(tabla_variograma, row.names = FALSE)

cat(sprintf("\n================ CORTE %d, SEMANA %02d ================\n", CORTE[1], CORTE[2]))
fila <- tabla_cortes[tabla_cortes$anio == CORTE[1] & tabla_cortes$semana == CORTE[2], ]
if (!nrow(fila)) {
  cat("No existe este corte. Cambie CORTE.\n")
} else {
  cat("Modelo ganador y resultado final:\n"); print(fila, row.names = FALSE)
  cat("\nPrimeros 12 candidatos validos, ordenados por RMSE LOO:\n")
  print(head(tabla_competencia[tabla_competencia$estado_candidato == "VALIDO",
                              c("modelo", "escala", "vif_max", "MAE", "RMSE", "R2")], 12), row.names = FALSE)
  cat("\nCandidatos excluidos por colinealidad:", sum(tabla_competencia$colineal), "de", nrow(tabla_competencia), "\n")
}
cat("\nTablas guardadas en:\n", SALIDA, "\n", sep = "")
