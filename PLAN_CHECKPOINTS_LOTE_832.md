# Plan de ejecución con checkpoints

El lote contiene 832 cortes: 16 años × 52 semanas ISO. Cada corte se procesa de forma independiente; un fallo en una semana no borra ni bloquea los resultados de las demás.

## C0 - Prueba de equivalencia

Antes del lote, ejecutar tres cortes ya conocidos: seco, intermedio y húmedo. El objetivo no es escoger modelos, sino comprobar que la versión final en R base reproduce mapas, semivariograma, MCO, MCG y predicción sin depender de `gstat` ni `spdep`.

## C1 - Primer año completo

Ejecutar 2010 completo: 52 cortes. Guardar resultados parciales cada 10 cortes y generar un reporte anual. Se revisa que no haya errores sistemáticos de covarianza, semivariograma o validación.

## C2 - Lote completo

Si C1 es estable, ejecutar 2011-2025. Se escribe un checkpoint cada 10 cortes y un resumen al cerrar cada año. No se generan 832 juegos de figuras: se guardan tablas para todos y figuras detalladas para C0, un corte ganador representativo por año y los cortes marcados para revisión.

## Reglas automáticas por corte

- Control de calidad fallido: `FALLIDO_CONTROL`, guardar motivo y continuar.
- VIF alto: excluir solo ese modelo candidato.
- Nugget puro o Moran residual no significativo: finalizar MCO.
- Semivariograma no ajustable: `REVISAR_VARIAGRAMA`, conservar MCO y continuar.
- Matriz `Sigma` no positiva definida: `REVISAR_COVARIANZA`, no aplicar MCG/kriging.
- Kriging no mejora validación: conservar MCO/MCG.
- Kriging mejora: guardar predicción y varianza de predicción.

## Reglas automáticas de pausa del lote

La ejecución se detiene para revisión técnica si ocurre cualquiera de estas condiciones dentro de un año:

1. Más del 5% de cortes falla control de calidad.
2. Más del 10% termina con error de código, semivariograma no ajustable o `Sigma` no válida.
3. Aparece una predicción no finita, negativa extrema o con varianza negativa no atribuible a redondeo.

Que gane MCO, MCG o kriging no pausa el lote: son resultados estadísticos válidos. Solo se pausa ante señales de un problema de datos o implementación.

## Archivos que se guardan

- `resumen_decisiones.csv`: una fila por corte y decisión final.
- `modelos_candidatos.csv`: métricas y VIF de los 57 candidatos por corte.
- `semivariogramas.csv`: parámetros y ajuste MCO/MCP de exponencial, esférico y gaussiano.
- `mco_mcg.csv`: coeficientes y errores estándar de ambas estimaciones.
- `validacion.csv`: MAE, RMSE, sesgo, R² predictivo y cobertura.
- `predicciones.rds`: predicción y varianza por celda solo para el modelo seleccionado.
- `figuras/`: únicamente cortes C0, representativos anuales y alertas.
