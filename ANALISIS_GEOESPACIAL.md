# Análisis Geoespacial — Precipitación en el Valle del Cauca

Este documento explica y justifica el análisis geoespacial realizado en
`02_Analisis_Geoespacial.R`, que extiende el trabajo de reconstrucción de
datos y EDA descriptivo (`01_Preparar_Verificar_Dataset.R`, `EDA_final.R`)
hacia el análisis propiamente geoestadístico. El dataset tiene **dos
dimensiones que se analizan por separado y luego se conectan**: una
dimensión **espacial** (686 celdas que cubren el Valle del Cauca) y una
dimensión **temporal** (16 años × 52 semanas, 2010–2025).

---

## 1. Por qué se necesita este análisis (y no basta con el EDA)

El EDA descriptivo (`EDA_final.R`) ya mostró que precipitación, temperatura,
radiación y altitud están espacialmente correlacionadas entre sí (matriz de
correlación, mapas promedio por celda) y que la relación con precipitación
cambia según la escala de agregación (pooled, espacial, temporal, anomalías).
Sin embargo, ese análisis **no cuantifica formalmente la dependencia
espacial** ni verifica si esa dependencia sigue presente después de ajustar
una tendencia — dos preguntas centrales de la teoría de variables
regionalizadas (Matheron) que sustenta el uso de kriging:

1. **¿La precipitación de una celda depende de sus vecinas más de lo que
   dependería si estuviera distribuida al azar en el espacio?** (autocorrelación
   espacial — Ley de Tobler: "todo está relacionado con todo, pero las cosas
   cercanas están más relacionadas que las distantes").
2. **Si se remueve una tendencia de gran escala (gradiente espacial +
   altitud), ¿queda estructura espacial en lo que sobra (el residual)?**
   Si la respuesta es sí, un modelo puramente determinístico (solo la
   tendencia) es insuficiente, y se necesita un componente estocástico
   espacial (semivariograma + kriging) para capturar esa estructura fina.

Estas dos preguntas se responden en las secciones 1 y 2 del script, **antes**
de ajustar cualquier semivariograma — siguiendo el mismo orden lógico que
plantean los ejemplos de geoestadística de clase.

---

## 2. Autocorrelación espacial: Índice de Moran e Índice de Geary

**Método:** se usa exactamente la metodología de clase (`funciones_auxiliares.R`,
sin `spdep`/`gstat`): una matriz de pesos de vecindad tipo "reina" (queen
contiguity) construida con `terra::adjacent()`, fila-estandarizada, y las
fórmulas cerradas de Moran (`I`) y Geary (`C`) por álgebra matricial.

**Por qué "reina" y no "torre":** en una grilla regular como esta (celdas de
~5.57 km), la contigüidad reina incluye también las 4 celdas diagonales, lo
que es más apropiado para un fenómeno continuo como la lluvia (que no tiene
una razón física para respetar solo las direcciones cardinales).

**Justificación del muestreo temporal (4 semanas):** en vez de calcular Moran
solo para una semana (lo cual respondería la pregunta espacial pero no la
temporal), se calculó para **4 semanas representativas de las dos temporadas
climáticas del régimen bimodal del Valle del Cauca** (lluvia: mayo y
noviembre; seca: febrero y agosto), en años distintos, para verificar si la
autocorrelación espacial es una propiedad estable del fenómeno o depende de
la época del año.

**Resultado:**

| Año-semana | Temporada | Moran I | Geary C |
|---|---|---:|---:|
| 2015-s18 | lluvia (may) | 0.950 | 0.040 |
| 2020-s44 | lluvia (nov) | 0.968 | 0.035 |
| 2016-s05 | seca (feb) | 0.762 | 0.245 |
| 2021-s32 | seca (ago) | 0.979 | 0.025 |

El valor esperado de Moran I bajo ausencia de autocorrelación es
`-1/(n-1) ≈ -0.0015`. Las 4 semanas están muy por encima de ese valor
(0.76–0.98), con Geary C cercano a 0 (máxima autocorrelación positiva)
— **la autocorrelación espacial es fuerte y está presente en las dos
temporadas climáticas**, no es un artefacto de una sola época del año.
La semana más débil (2016-s05, seca) también es la de menor precipitación
media (2.6 mm) — con casi toda el área seca, hay menos variación real que
autocorrelacionar, lo cual es coherente, no contradictorio.

---

## 3. Tendencia espacial y autocorrelación del residual

**Método:** para cada una de las 4 semanas se ajusta la misma tendencia
lineal que ya se usaba en la fase de modelación (`02_Modelacion.R`):

```
precipitacion_mm ~ x_km + y_km + altitud_m
```

(coordenadas convertidas a kilómetros con una proyección local simple,
apropiada para el tamaño del Valle del Cauca — no se necesita una
proyección cartográfica completa para un área de este tamaño). Después se
calcula el Índice de Moran **sobre los residuales** de esa tendencia, no
sobre el dato crudo.

**Por qué este paso es indispensable:** si el residual no tuviera
autocorrelación espacial, la tendencia ya habría capturado toda la
estructura espacial del fenómeno, y ajustar un semivariograma sería
innecesario (bastaría con la regresión). Comprobar lo contrario es lo que
justifica metodológicamente el uso de kriging en este proyecto.

**Resultado:**

| Año-semana | R² tendencia | Moran I crudo | Moran I residual |
|---|---:|---:|---:|
| 2015-s18 | 0.721 | 0.950 | 0.855 |
| 2020-s44 | 0.663 | 0.968 | 0.897 |
| 2016-s05 | 0.035 | 0.762 | 0.752 |
| 2021-s32 | 0.794 | 0.979 | 0.893 |

La tendencia reduce el Índice de Moran en todas las semanas, pero **nunca
lo lleva cerca de cero** — el residual sigue teniendo autocorrelación
espacial fuerte (0.75–0.90). Esto confirma que un modelo de tendencia
determinístico por sí solo es insuficiente, y justifica formalmente pasar
al semivariograma y kriging (sección 4). El caso más extremo es 2016-s05
(seca), donde el R² de la tendencia es casi nulo (0.035) — en una semana
casi sin lluvia en toda el área, ni la geografía ni la altitud explican
mucho, porque casi todo el departamento comparte el mismo valor (cerca de
cero); ahí el semivariograma y el kriging cargan con casi todo el trabajo
de la predicción.

---

## 4. Semivariograma y kriging universal, por semana

**Método:** idéntico al de `02_Modelacion.R` (para no introducir una
metodología nueva y no validada): semivariograma empírico por bins de
distancia (12 bins, corte en 180 km), modelo exponencial ajustado por
mínimos cuadrados ponderados por el número de pares en cada bin, kriging
universal resolviendo el sistema lineal aumentado con el multiplicador de
Lagrange, y validación cruzada leave-one-out (LOO) sobre una muestra de 120
celdas por semana (para mantener el tiempo de cómputo razonable sin perder
representatividad).

**Resultado:**

| Año-semana | Temporada | Nugget | Sill | Rango (km) | RMSE (mm) | RMSE relativo |
|---|---|---:|---:|---:|---:|---:|
| 2015-s18 | lluvia (may) | 0 | 240.9 | 142.2 | 2.61 | 6.3% |
| 2020-s44 | lluvia (nov) | 0 | 954.9 | 92.0 | 5.15 | 7.8% |
| 2016-s05 | seca (feb) | 0 | 95.1 | 27.5 | 3.59 | **137.9%** |
| 2021-s32 | seca (ago) | 0 | 536.8 | 76.8 | 5.96 | 16.6% |

**Interpretación y justificación de por qué se reporta así:**

- **El nugget sale en 0 en las 4 semanas.** Esto es razonable para datos de
  satélite/reanálisis interpolados a una grilla común (CHIRPS y NASA POWER
  reconstruido): a diferencia de estaciones puntuales con error de medición
  instrumental, esta grilla no tiene "error de pepita" por duplicados en el
  mismo sitio — el semivariograma parte casi de cero en distancia cero.
- **El sill y el rango cambian mucho de una semana a otra** (sill de 95 a
  955; rango de 27 a 142 km) — **la estructura espacial de la precipitación
  no es fija en el tiempo**, depende del evento climático puntual de esa
  semana. Esto es un hallazgo relevante para el proyecto: un solo
  semivariograma "promedio" no representaría bien ninguna semana en
  particular.
- **El RMSE relativo de 137.9% en la semana seca (2016-s05) no es un error
  del modelo, es una limitación esperable:** con precipitación media de solo
  2.6 mm esa semana, errores absolutos pequeños (3.59 mm) se convierten en
  porcentajes de error enormes. Se reporta así, sin ocultarlo, porque es
  información real sobre las limitaciones del modelo en semanas casi secas
  — coherente con lo ya documentado en el EDA sobre el 9.9% de semanas con
  precipitación cero.

---

## 5. Autocorrelación temporal (dimensión tiempo)

**Método:** hasta este punto todo el análisis trata cada semana como una
superficie espacial independiente. Para revisar explícitamente la dimensión
**temporal** (si el valor de una celda en una semana predice su propio valor
en la semana siguiente), se calculó la anomalía de precipitación
(`precipitacion_mm - clim_precipitacion_mm`, quitando el ciclo estacional) y
su función de autocorrelación (ACF) hasta 12 semanas de rezago, en una
muestra aleatoria de 6 celdas (semilla fija para reproducibilidad).

**Por qué se usa la anomalía y no el dato crudo:** el dato crudo tiene
autocorrelación temporal trivial por el ciclo estacional (las semanas de
temporada de lluvia se parecen entre sí simplemente porque ambas son
temporada de lluvia, no porque una cause a la otra). Restar la climatología
aísla la variación semana a semana que no se explica por la estación del
año — la pregunta relevante para decidir si hay "memoria" temporal real.

**Resultado:**

| Celda | ACF rezago 1 | ACF rezago 2 |
|---|---:|---:|
| 633 | 0.067 | 0.003 |
| 294 | 0.024 | 0.039 |
| 557 | **0.326** | **0.275** |
| 623 | 0.120 | 0.052 |
| 108 | 0.197 | 0.028 |
| 164 | 0.050 | -0.004 |

En 5 de las 6 celdas, la autocorrelación a una semana de rezago es baja
(0.02–0.20) y cae rápido — la anomalía de precipitación se comporta cerca
de un proceso sin memoria de una semana a la siguiente. **La celda 557 es
una excepción clara**, con autocorrelación sostenida alrededor de 0.3 en
varios rezagos (no decae hacia cero) — sugiere que en esa ubicación
específica los episodios húmedos o secos tienden a persistir más de una
semana. No se generaliza esto a todo el departamento con solo 6 celdas de
muestra; se deja documentado como un patrón a explorar con más celdas si el
proyecto lo requiere.

**Justificación de la decisión metodológica que esto habilita:** dado que la
autocorrelación temporal es en general baja, se justifica el enfoque
adoptado en las secciones 2–4 (modelar cada semana como una realización
espacial separada, sin un componente espacio-temporal conjunto explícito) —
no porque el tiempo sea irrelevante (la sección 4 ya mostró que los
parámetros del semivariograma cambian con la época del año), sino porque una
semana no predice fuertemente a la siguiente, así que no hay una ganancia
clara de modelar la dependencia temporal semana a semana además de la
espacial.

---

## 6. Limitaciones y decisiones pendientes

- El muestreo de 4 semanas (2 lluvia, 2 seca) es representativo pero no
  exhaustivo; los patrones de sill/rango podrían variar más si se
  incluyeran más semanas o los años extremos ya identificados en el EDA
  (2010–2013, época de La Niña).
- La ACF temporal se calculó sobre solo 6 celdas de muestra; el hallazgo de
  la celda 557 debe tratarse como una observación puntual, no como una
  conclusión general del departamento.
- El nugget en 0 debe revisarse si en el futuro se trabaja con estaciones
  puntuales reales (no interpoladas), donde sí es esperable observar un
  nugget positivo por variabilidad de micro-escala y error de instrumento.
- Sigue pendiente decidir formalmente entre las variantes vecino/bilineal de
  radiación (ver `README.md` del proyecto) — este análisis geoespacial se
  hizo con las variantes bilineales por consistencia con el resto del
  proyecto, pero no zanja esa decisión.

## 7. Archivos generados

Todos los resultados numéricos de este análisis se guardan en
`analisis_geoespacial_resultados/`:

- `moran_geary_por_semana.csv`
- `moran_crudo_vs_residual.csv`
- `resumen_semivariograma_kriging.csv`
- `autocorrelacion_temporal_acf.csv`

Los gráficos (panel de semivariogramas, panel de ACF) se muestran en la
pestaña Plots de RStudio al correr `02_Analisis_Geoespacial.R`, sin guardado
automático a PDF/PNG, siguiendo el mismo criterio ya usado en `EDA_final.R`.
