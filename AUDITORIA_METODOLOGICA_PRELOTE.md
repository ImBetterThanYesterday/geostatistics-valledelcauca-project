# Auditoría metodológica antes del lote completo

Esta revisión contrasta el proyecto, los ejemplos R y los materiales de clase agregados el 12 de septiembre de 2026. La unidad de análisis definitiva es una superficie de precipitación de un año y una semana ISO, con 686 píxeles. No se promedian semanas ni años.

## Fases y estado

| Fase | Estado | Regla final |
|---|---|---|
| 1. Control de calidad | Listo | Exigir 686 celdas únicas, sin NA, precipitación no negativa y covariables completas. |
| 2. Soporte espacial | Listo | Plantilla CHIRPS de 0.05°, píxel cuyo centro está dentro del Valle; coordenadas convertidas a km antes de distancias. |
| 3. EDA | Listo, se ejecuta por corte | Posting/mapa, histograma, dispersión contra coordenadas y covariables, distancias y extremos. |
| 4. Media | Revisado | Comparar la escalera completa de modelos ambientales y de tendencia. Aplicar VIF a modelos con más de un predictor. |
| 5. MCO y residuales | Pendiente de versión final | Elegir por error LOO en mm y calcular residuales del candidato ganador. |
| 6. Dependencia residual | Pendiente de versión final | Moran I manual y lectura del semivariograma. |
| 7. Semivariograma | Pendiente de versión final | Construir sobre residuales, con lags documentados, número de pares por lag y diagnóstico direccional. |
| 8. Ajuste teórico | Pendiente de versión final | Exponencial, esférico y gaussiano por MCO y MCP con múltiples valores iniciales; guardar estabilidad de parámetros. |
| 9. Covarianza y MCG | Incorporado al algoritmo final | Construir Sigma desde el variograma elegido y reestimar el mismo modelo de media por MCG/GLS; comparar con MCO. |
| 10. Predicción | Pendiente de versión final | M0 lleva a kriging ordinario; cualquier modelo con tendencia lleva a kriging universal manual, usando la media MCG. |
| 11. Validación | Pendiente de versión final | LOO y RMSE, con predicciones en mm. |
| 12. Comparación final | Pendiente de versión final | Comparar MCO contra MCG/kriging y conservar la alternativa que mejore de forma verificable. |

## Modelos de media que entran al lote

| Nombre | Modelo |
|---|---|
| M0 | `precipitacion ~ 1` |
| Mxy | `precipitacion ~ x_km + y_km` |
| Malt | `precipitacion ~ altitud_m` |
| MT | `precipitacion ~ temperatura_bilineal_c` |
| MR | `precipitacion ~ radiacion_bilineal_mj_m2_dia` |
| Mxy_alt | `precipitacion ~ x_km + y_km + altitud_m` |

### Combinaciones ambientales

| Nombre | Modelo |
|---|---|
| MTR | `precipitacion ~ temperatura + radiacion` |
| MTA | `precipitacion ~ temperatura + altitud` |
| MRA | `precipitacion ~ radiacion + altitud` |
| MTRA | `precipitacion ~ temperatura + radiacion + altitud` |

### Tendencia espacial lineal más ambiente

| Nombre | Modelo |
|---|---|
| Mxy_T | `precipitacion ~ x + y + temperatura` |
| Mxy_R | `precipitacion ~ x + y + radiacion` |
| Mxy_TR | `precipitacion ~ x + y + temperatura + radiacion` |
| Mxy_TA | `precipitacion ~ x + y + temperatura + altitud` |
| Mxy_RA | `precipitacion ~ x + y + radiacion + altitud` |
| Mxy_TRA | `precipitacion ~ x + y + temperatura + radiacion + altitud` |

La lista contiene 16 modelos de media. Cada uno se ajusta con respuesta original, raíz cuadrada y `log1p`: 48 candidatos por corte. VIF > 10 descarta solo ese candidato, no el corte completo. El número de condición y la validación LOO complementan la decisión.

No se usan términos cuadráticos `x²`, `y²` ni `x·y`: no se implementaron como parte de los ejemplos aplicados del profesor y el proyecto exige utilizar únicamente metodologías trabajadas en clase.

## Semivariogramas que entran al lote

- Exponencial.
- Esférico.
- Gaussiano.

Matérn aparece en las diapositivas teóricas, pero no se incluye en la selección automática: requiere un parámetro adicional de suavidad y no está en el ejercicio aplicado. La anisotropía se revisa con variogramas direccionales y exige revisión explícita antes de ajustar un modelo anisotrópico.

## Controles añadidos desde el material de clase

- Distancias: se calculan sobre `x_km`, `y_km`, nunca directamente sobre grados de longitud y latitud.
- Lags: cada fila del semivariograma guarda distancia media y número de pares. Lags con pocos pares se señalan para revisión.
- Ajuste: cada familia se intenta desde varios rangos iniciales derivados de cuantiles de las distancias. Se comparan MCO y MCP para el ajuste de parámetros del semivariograma.
- Nugget puro: si la contribución espacial queda prácticamente en cero y no hay mejora frente al modelo nugget, el corte termina en MCO y no se fuerza kriging.
- Covarianza: antes de MCG se comprueba que `Sigma` sea simétrica y positiva definida. Si falla, se registra `REVISAR_COVARIANZA` y no se genera una predicción engañosa.
- Incertidumbre: para los modelos que terminan en kriging se guarda también la varianza de predicción, no solo el mapa predicho.

## Paso MCO → MCG/GLS

1. MCO estima la tendencia y deja residuales.
2. El semivariograma de los residuales define nugget, sill y rango.
3. Esos parámetros construyen la matriz de covarianza `Sigma` entre las 686 celdas.
4. MCG/GLS vuelve a estimar **la misma fórmula ganadora**, pero reconoce que los errores cercanos se parecen entre sí.
5. Se guardan coeficientes, errores estándar y ajuste de MCO y MCG. El kriging incorpora posteriormente la corrección local de los residuales.

MCG no agrega covariables ni cambia de M0, MT o Mxy_TRA: cambia la forma de estimar sus coeficientes cuando ya sabemos que los errores no son independientes.

## Evidencia de clase usada

- `Actividad_2.pdf`, sección de precipitación: exige comparar una semana sin promediar, tendencia MCO, variograma de residuales, ajuste exponencial con inicios múltiples, comparación exponencial/esférico/gaussiano, kriging universal manual y LOO.
- `Actividad_1.pdf`: muestra una respuesta ambiental en función de dos covariables y el paso de MCO a MCG/GLS cuando el error tiene dependencia espacial.
- `Ejemplo3_Geoestadística.R`: usa `precip_media ~ x_km + y_km`, residuales, modelo exponencial y kriging universal.
- `Practic_Geostatistics.R`: usa `Z ~ 1` y kriging ordinario.

## Restricción de implementación

El enunciado del proyecto permite librerías adicionales solo para gráficos y carga de mapas/Excel. Por tanto, `gstat` y `spdep`, usados en el checkpoint exploratorio, no pueden formar parte del archivo final de entrega. La implementación final debe calcular con R base la matriz de distancias, Moran, semivariograma, ajuste por `optim`, matriz de covarianza, MCG y kriging, igual que los ejemplos del profesor.
