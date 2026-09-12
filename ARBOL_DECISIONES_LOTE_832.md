# Árbol de decisiones: 832 modelos espaciales de precipitación

Unidad del análisis: un corte de **un año y una banda semanal**, con 686 celdas. El lote no promedia años ni bandas y no usa climatologías. La respuesta única siempre es `precipitacion_mm`; temperatura, radiación, altitud y coordenadas solo pueden entrar como covariables de tendencia.

```text
INICIO: año + banda
  |
  |-- Control de calidad: 686 celdas únicas, sin NA, precipitación >= 0
  |     |-- No cumple -> FALLIDO_CONTROL; guardar motivo; siguiente corte
  |     '-- Cumple
  |
  |-- Ajustar modelos de media definidos para el proyecto y tres escalas de respuesta
  |     |-- VIF > 10 -> ese candidato se excluye por colinealidad
  |     '-- Elegir el candidato con menor RMSE LOO en mm
  |
  |-- Moran sobre residuales del modelo ganador
  |     |-- p >= 0.05 -> FINAL_MCO
  |     '-- p < 0.05 -> dependencia espacial residual
  |
  |-- Diagnóstico direccional (alerta de anisotropía; no fuerza ajuste aún)
  |
  |-- Semivariograma residual isotrópico
  |     |-- insuficiente/no converge -> REVISAR_VARIAGRAMA
  |     '-- Ajustar exponencial, esférico y gaussiano
  |             '-- Elegir menor suma de cuadrados con parámetros admisibles
  |
  |-- Construir matriz de covarianza y reestimar por MCG/GLS
  |     |-- Comparar MCO vs MCG: coeficientes, errores estándar y ajuste
  |     '-- Usar la tendencia MCG como componente de media espacial
  |
  |-- Kriging LOO con variograma fijo
  |     |-- predicciones no finitas o extremas -> FINAL_MCO_KRIGING_INESTABLE
  |     |-- RMSE no mejora al MCO al menos 1% -> FINAL_MCO
  |     '-- mejora >= 1%
  |             |-- M0 ganador -> FINAL_KRIGING (ordinario)
  |             '-- otro modelo ganador -> FINAL_KRIGING (universal)
  |
  '-- Guardar: variables, VIF, métricas MCO/kriging, Moran, alerta direccional,
      variogramas candidatos, decisión y motivo.
```

## Modelos de media candidatos aprobados para el lote

| Grupo | Modelos |
|---|---|
| Referencia y un predictor | M0, Mxy, Malt, MT, MR |
| Ambiente combinado | MTR, MTA, MRA, MTRA |
| Tendencia lineal + ambiente | Mxy_alt, Mxy_T, Mxy_R, Mxy_TR, Mxy_TA, Mxy_RA, Mxy_TRA |

Las fórmulas exactas están en `AUDITORIA_METODOLOGICA_PRELOTE.md`. No se incluyen tendencias cuadráticas porque no aparecen en los ejemplos aplicados del curso. La comparación usa la misma métrica para todas las escalas: predicción devuelta a milímetros y evaluación LOO. La transformación no es una covariable ni una “pista” de precipitación: es otra escala en la que se ajusta la respuesta. VIF > 10 descarta solamente la combinación colineal, no el corte.

## Semivariogramas aprobados para selección automática

Se comparan el **exponencial**, **esférico** y **gaussiano**, tal como exige la Actividad 2. Se elige el de menor suma de cuadrados entre el semivariograma empírico y teórico, siempre que nugget, sill parcial y rango sean admisibles. Matérn se presenta en las diapositivas como una familia teórica, pero no se automatiza en el lote: requiere elegir además un parámetro de suavidad y no aparece en el ejercicio aplicado. La anisotropía se diagnostica con variogramas direccionales; solo se ajusta tras una revisión específica, no por una regla ciega.

## Límites iniciales, sujetos al checkpoint

- VIF máximo permitido: 10.
- Moran: p < 0.05 para intentar el componente espacial.
- Mejora mínima de kriging: 1% de RMSE respecto al MCO ganador.
- Predicción inestable: algún valor no finito o una predicción superior a tres veces el máximo observado del mismo corte. Se conserva el MCO y se deja trazabilidad para revisar el caso.
- Anisotropía: alerta exploratoria cuando el cociente direccional de semivarianza de corto alcance supera 2. No se ajusta automáticamente hasta validar que esa señal sea estable.

## Nota de validez

La validación de kriging de este primer lote es LOO con variograma fijo, por lo que es útil como filtro y comparación interna, pero es optimista. Antes de declarar resultados finales se deben repetir los modelos seleccionados con una validación más estricta que vuelva a estimar el variograma dentro de cada partición. La versión de entrega debe calcular distancias, semivariogramas, covarianza, MCG y kriging manualmente con R base, siguiendo las condiciones del proyecto.
