# Bitácora metodológica — Proyecto 1

**Proyecto:** precipitación semanal en el Valle del Cauca  
**Curso:** Analítica de Datos, 2026-2S  
**Estado:** reinicio desde la auditoría  
**Última actualización:** 2026-09-04

## Propósito

Esta es la fuente oficial para registrar hechos verificados, decisiones
aprobadas, alternativas, riesgos, evidencia reproducible y preguntas
pendientes. Una observación no se convierte automáticamente en decisión y una
hipótesis no se presenta como hecho hasta verificarla.

### Convenciones

| Código | Significado |
|---|---|
| `H-###` | Hecho verificado con datos o código |
| `D-###` | Decisión aprobada |
| `P-###` | Propuesta pendiente o alternativa |
| `R-###` | Riesgo o limitación |
| `A-###` | Acción siguiente |

## Reglas de trabajo

1. Los archivos originales nunca se sobrescriben.
2. Toda limpieza conserva el valor y la categoría originales.
3. Los resultados se muestran en consola y en **Plots**, y se guardan en disco.
4. Antes de imputar se distingue el tipo de ausencia.
5. Antes de modelar se fija el soporte espacial y temporal.
6. Las decisiones importantes se acuerdan antes de modificar los datos.
7. El rendimiento se evalúa fuera de muestra, preferiblemente por bloques.
8. Una decisión requiere evidencia, no solo apariencia visual.
9. Durante el desarrollo usamos módulos; la entrega será un único R autónomo.
10. El análisis principal respetará los procedimientos vistos en clase.

---

# 1. Inventario y linaje

| Variable | Fuente | Producto anual | Climatología |
|---|---|---|---|
| Precipitación | CHIRPS | 52 bandas/año | 52 bandas |
| Temperatura a 2 m | NASA POWER | Primeras 52 bandas POWER | 52 bandas |
| Radiación solar | NASA POWER | Últimas 52 bandas POWER | 52 bandas |
| Altitud | SRTM | Una capa estática | No aplica |

Periodo: 2010–2025, nominalmente 832 semanas por variable dinámica.

Una climatología semanal es el promedio de una misma semana entre años:

\[
\bar Y_s(u)=\frac{1}{16}\sum_{a=2010}^{2025}Y_{a,s}(u).
\]

No es una medición independiente ni una segunda fuente.

---

# 2. Hechos verificados

## H-001 — Geometría

Los rásteres inspeccionados tienen 39 filas, 37 columnas, 1,443 celdas,
resolución aproximada de 0.05° y CRS WGS84. Compartir geometría no implica
compartir cobertura válida.

## H-002 — Reglas de máscara

- Centro del píxel dentro de GADM: 686 celdas.
- Cualquier intersección con GADM: 789 celdas.
- Diferencia: 103 celdas fronterizas.

En esas 103 celdas, la fracción mediana dentro del Valle es 19%; 41 tienen
menos de 10% dentro y solo 5 tienen al menos 50% dentro.

## H-003 — Cobertura interior semanal

Sobre la máscara de 686 celdas:

| Variable | Válidos | NA internos | Cobertura |
|---|---:|---:|---:|
| Precipitación | 686 | 0 | 100% |
| Altitud | 686 | 0 | 100% |
| Temperatura | 122 | 564 | 17.78% |
| Radiación | 27 | 659 | 3.94% |

Los conteos son constantes en 2010–2025. La ausencia POWER es espacial y
sistemática, no una pérdida temporal ocasional.

## H-004 — Intersección cruda

En 2024-S10, solo una celda interior contiene simultáneamente precipitación,
temperatura, radiación y altitud sin imputar. No permite estimar un modelo.

## H-005 — Códigos CHIRPS

Dentro de la máscara por centro no hay valores negativos. Códigos extremos como
`-69993` y `-89991` aparecen en celdas fronterizas de la máscara por
intersección y no representan precipitación física.

## H-006 — Altitud negativa pequeña

Existe un valor SRTM cercano a −2.9 m. No se declaró inválido: puede ser posible
en costa o reflejar incertidumbre vertical. Se conserva con observación.

## H-007 — Climatologías derivadas

Las climatologías se reprodujeron promediando los 16 valores anuales de cada
semana:

| Variable | Correlación | Diferencia máxima absoluta |
|---|---:|---:|
| Temperatura | 1.000000 | 0.00000095 |
| Radiación | 1.000000 | 0.00000095 |

Por tanto, heredan exactamente los NA de los archivos anuales.

## H-008 — Franjas POWER

Temperatura y radiación aparecen en franjas, no como superficies completas.
Esto contradice la descripción del enunciado de capas remuestreadas comparables
píxel a píxel.

**Hipótesis pendiente:** pudieron insertarse nodos POWER en una plantilla CHIRPS
sin completar correctamente el remuestreo. Aún no se presenta como causa
demostrada.

## H-009 — Receta semanal de radiación casi reproducida

**Fecha:** 2026-09-04

En dos nodos independientes se consultó NASA POWER diario para todo 2024 con:

- parámetro `ALLSKY_SFC_SW_DWN`;
- comunidad `AG`;
- tiempo `LST`;
- unidad `MJ/m²/día`;
- bloques de siete días desde el 1 de enero;
- promedio de los siete valores diarios.

El resultado coincide con el GeoTIFF suministrado en 51 de 52 semanas en ambos
nodos, con correlaciones 0.9975 y 0.9952. Para 2024-S10 las coincidencias son
exactas: 18.31571 y 14.83714 MJ/m²/día.

La excepción es la semana 1, con diferencia cercana a 1 MJ/m²/día. Puede
corresponder a una revisión posterior del archivo fuente o a una particularidad
de preparación del primer intervalo. Todavía debe investigarse.

## H-010 — Temperatura reproduce el patrón, no los valores exactos

Con `T2M`, comunidad `AG`, `LST` y promedio en los mismos bloques semanales:

| Nodo | Correlación | MAE | RMSE | Sesgo entregado − NASA |
|---|---:|---:|---:|---:|
| punto 1 | 0.9779 | 0.3950 °C | 0.4085 °C | −0.3950 °C |
| punto 2 | 0.9918 | 0.1027 °C | 0.1478 °C | −0.0968 °C |

La receta temporal parece correcta, pero la diferencia cambia espacialmente.
Debe comprobarse si proviene de una versión histórica de NASA, de la celda
nativa asignada o del remuestreo usado al crear el GeoTIFF.

## R-001 — Versionado de la fuente NASA

NASA puede revisar productos históricos. Una reconstrucción actual puede seguir
el mismo patrón y metodología sin ser idéntica a un archivo descargado en otra
fecha. Por eso se conservan los CSV crudos y los metadatos de cada consulta.

---

# 3. Clasificación permanente de calidad

| Categoría | Definición | Tratamiento |
|---|---|---|
| `Valido` | Dato admisible dentro del soporte | Conservar |
| `NA_estructural` | Ausencia fuera del Valle | Nunca imputar |
| `NA_interno_original` | NA que ya existía dentro | Investigar |
| `NA_control_calidad` | Valor invalidado por una regla documentada | Reparar solo si se justifica |
| `Dato_fuera_Valle` | Valor almacenado fuera del soporte | Excluir del análisis principal |

Cada corrección futura conservará: valor original, categoría original, valor
corregido, método, indicador de imputación y versión del procesamiento.

---

# 4. Decisiones aprobadas

## D-001 — Máscara oficial por centro

**Fecha:** 2026-09-04 — **Estado:** aprobada

El soporte espacial objetivo contiene las 686 celdas cuyo centro está dentro
del polígono GADM del Valle.

**Justificación:** es reproducible, coherente con coordenadas puntuales, evita
celdas que apenas tocan el límite y coincide con CHIRPS-altitud completos.

## D-002 — No imputar POWER en fases 1–2

**Fecha:** 2026-09-04 — **Estado:** aprobada

Los huecos permanecen como `NA_interno_original`. Se retira del flujo reiniciado
la imputación IDW exploratoria. Completar 82.2% de temperatura y 96.1% de
radiación sin validar sería principalmente extrapolación.

## D-003 — Soporte maestro de 686 celdas

**Fecha:** 2026-09-04 — **Estado:** aprobada

El soporte objetivo no se reduce por la disponibilidad de POWER. Los NA se
registran sobre las 686 celdas y no se trabaja con la única celda completa.

---

# 5. Alternativas para POWER

## P-001 — Reconstrucción desde NASA

**Estado:** alternativa preferida; pendiente de prueba piloto.

```text
NASA POWER diario en resolución original
                  ↓
identificar unidades, tiempo y agregación
                  ↓
temperatura y radiación semanales
                  ↓
remuestreo a la plantilla CHIRPS
                  ↓
máscara aprobada de 686 celdas
```

Antes de descargar 2010–2025 se intentará reproducir valores válidos entregados.
Se probarán parámetro, comunidad `AG/RE`, unidad, LST/UTC, definición de semana,
promedio/suma y método de remuestreo.

**Criterios de aprobación:** coincidencia en varios lugares y fechas, cobertura,
unidades físicas, continuidad espacial y reproducibilidad.

## P-002 — Modelo base sin POWER

**Estado:** aceptado conceptualmente; pendiente en el flujo reiniciado.

Usará precipitación, altitud y componentes espaciales/temporales. Será el punto
de referencia para medir el aporte incremental de POWER reconstruido.

## P-003 — Interpolación de POWER incompleto

**Estado:** solo sensibilidad.

Podría compararse IDW, regresión, bilineal o kriging, pero la distribución en
franjas produce extrapolación extensa. No será el producto principal sin una
validación satisfactoria.

---

# 6. Estado de fases

| Fase | Estado | Producto |
|---|---|---|
| 1. Auditoría y calidad | Completada, abierta | `00_Auditoria_Calidad.R` |
| 2. Soporte común | Completada | `01_Soporte_Comun.R` |
| 3. EDA correcto | En espera | Depende de POWER |
| 4–12 | No iniciadas en el reinicio | Modelado anterior no definitivo |

Los resultados anteriores de M3, Moran, variograma, MCG y kriging utilizaron
POWER completado por IDW. Son aprendizaje técnico, no resultados definitivos.

---

# 7. Mapa del código y resultados

## `00_Auditoria_Calidad.R`

Lee sin modificar, compara máscaras, clasifica ausencias, audita 832 semanas,
imprime en consola y guarda CSV, PDF y GeoTIFF.

## `01_Soporte_Comun.R`

Aplica D-001/D-003, conserva NA conforme a D-002, compara fuentes y demuestra
que la intersección completa es insuficiente.

```text
reinicio_desde_cero/
├── auditoria_resultados/
│   ├── conteos_categorias_semana.csv
│   ├── cobertura_2010_2025.csv
│   ├── fraccion_area_dentro_valle.csv
│   └── graficos_auditoria_nulos.pdf
└── soporte_resultados/
    ├── diferencias_datasets.csv
    ├── disponibilidad_soporte.csv
    ├── mascara_maestra_686.tif
    ├── soporte_maestro.rds
    └── visual_diferencias_y_soporte.pdf
```

---

# 8. Próximo punto de decisión

## A-001 — Piloto forense NASA POWER

1. Seleccionar nodos válidos suministrados.
2. Consultar POWER diario.
3. Probar configuraciones de unidades, fechas y agregación.
4. Intentar reproducir valores semanales entregados.
5. Registrar diferencias y gráficos.
6. Decidir si se aprueba la reconstrucción 2010–2025.

**Estado:** primera etapa completada; receta de radiación casi confirmada y
temperatura parcialmente confirmada.

El piloto muestra resultados en consola y **Plots**, preserva las respuestas
crudas y guarda tablas comparativas. El siguiente control propuesto es descargar
la región nativa únicamente para 2024, remuestrearla a la plantilla CHIRPS y
compararla contra todos los nodos suministrados antes de aprobar 2010–2025.

**Archivos:** `02_Piloto_Reconstruccion_POWER.R` y `power_piloto/`.

---

# 9. Investigación regional NASA POWER 2024

## H-011 — Resolución nativa y construcción semanal

**Fecha:** 2026-09-04 — **Tipo:** hecho reproducido

Los NetCDF regionales diarios actuales tienen resoluciones nativas diferentes:

| Variable | Parámetro | Resolución nativa | Unidad |
|---|---|---:|---|
| Temperatura a 2 m | `T2M` | 0.625° × 0.5° | °C |
| Radiación incidente | `ALLSKY_SFC_SW_DWN` | 1° × 1° | MJ/m²/día |

Se usaron comunidad `AG`, horario `LST` y promedios de bloques consecutivos de
7 días desde el 1 de enero. En 2024, las 52 semanas abarcan los días 1–364; los
días 365–366 no pertenecen a una semana completa del esquema entregado.

**Consecuencia:** la cuadrícula de 0.05° no agrega información observada; es un
remuestreo de una fuente mucho más gruesa. Esto debe declararse en el informe y
considerarse al interpretar coeficientes, significancia y validación espacial.

## H-012 — Método que reproduce los datos suministrados

**Fecha:** 2026-09-04 — **Tipo:** hallazgo

Comparación contra todas las celdas no nulas suministradas dentro del soporte:

| Variable | Método | Periodo | RMSE | Correlación |
|---|---|---|---:|---:|
| Temperatura | vecino | S02–S52 | 0.5052 °C | 0.9820 |
| Temperatura | bilineal | S02–S52 | **0.2446 °C** | **0.9959** |
| Temperatura | cúbico | S02–S52 | 0.7941 °C | 0.9801 |
| Radiación | vecino | S02–S52 | **0 exacto** | **1.0000** |
| Radiación | bilineal | S02–S52 | 0.1643 | 0.9937 |
| Radiación | cúbico | S02–S52 | 0.3984 | 0.9661 |

La coincidencia exacta de radiación en semanas 2–52 demuestra que la fuente
entregada fue construida con la misma agregación temporal y vecino más cercano.
La semana 1 presenta una discrepancia aislada cercana a 1 MJ/m²/día, compatible
con una revisión histórica del producto o una anomalía en la primera semana.

## H-013 — La diferencia de temperatura es espacial y estable

**Fecha:** 2026-09-04 — **Tipo:** diagnóstico

El candidato bilineal actual tiene RMSE 0.2559 °C usando S01–S52. Al restar solo
para diagnóstico el error medio fijo de cada celda, el RMSE cae a 0.0797 °C y el
MAE a 0.0295 °C. Quitar el promedio de cada semana apenas reduce el RMSE a
0.2441 °C. Además, excepto en S01, el error medio semanal oscila cerca de cero.

**Interpretación:** la temporalidad y la agregación semanal están bien
reconstruidas. La diferencia restante corresponde principalmente a un patrón
espacial persistente. Las explicaciones plausibles son una versión histórica de
NASA/MERRA-2 distinta a la disponible hoy o una alineación espacial anterior
ligeramente diferente. No hay evidencia para inventar valores ni para calibrar
con los pocos píxeles suministrados.

**Decisión cautelar:** no aplicar la corrección por celda. Haber demostrado que
reduce el error sirve para diagnosticar la causa; usarla sería sobreajustar la
reconstrucción a franjas incompletas y no tendría respaldo fuera de ellas.

## A-002 — Propuesta antes de descargar 2010–2025

**Estado:** pendiente de aprobación conjunta

1. Temperatura: bilineal como reconstrucción principal.
2. Radiación: vecino más cercano como producto fiel al proceso entregado.
3. Radiación bilineal: conservar como análisis de sensibilidad, porque es más
   suave, pero no representa valores nuevos observados.
4. Marcar S01-2024 como discrepancia de procedencia, sin reemplazarla mediante
   una corrección construida a mano.
5. Descargar 2010–2025 solo después de aprobar este criterio y automatizar la
   validación por año.

**Archivos:** `03_Reconstruccion_Regional_2024.R` y
`power_regional_2024/resultados/`.

## H-014 — Auditoría de procedencia dentro de `datos_proyecto_1`

**Fecha:** 2026-09-04 — **Tipo:** hallazgo de procedencia

Se inventariaron recursivamente los archivos sin modificarlos. La carpeta sí
contiene insumos diarios intermedios de CHIRPS: 192 GeoTIFF mensuales, uno por
mes entre 2010 y 2025, cuyas bandas son días. Estos archivos ya están recortados
y escritos en la cuadrícula de 0.05°, por lo cual son anteriores a los productos
semanales, aunque no necesariamente son las descargas crudas originales.

Para NASA POWER no se encontraron CSV, JSON, NetCDF ni GeoTIFF diarios nativos.
Solo están los 16 productos anuales `power_semanal_valle_<año>.tif`, cada uno con
104 bandas, y las dos climatologías derivadas. Los archivos POWER anuales fueron
creados consecutivamente como un lote el 2026-08-31 aproximadamente entre
20:12:24 y 20:13:02. Esto apoya que son salidas de un proceso común, pero no
permite recuperar el código, consulta, versión o archivo diario que los generó.

`datos_procesados/rasters_limpios.rds` tampoco contiene fuentes nativas: agrupa
los mismos productos semanales en pilas de 832 capas y registra
`modo_power = "combinado_104"`. El historial `.Rapp.history` encontrado solo
carga el límite GADM y no documenta la construcción de POWER.

**Conclusión:** `datos_proyecto_1` es la referencia oficial entregada y contiene
la cadena diaria útil de CHIRPS, pero no contiene los originales diarios de
POWER. Por tanto, para POWER debemos conservar los GeoTIFF como evidencia de
comparación y reconstruir desde NASA, salvo que el profesor pueda facilitar los
NetCDF/CSV originales o el script de preparación.

**Pregunta adicional para el profesor:** ¿existen los archivos diarios nativos
de NASA POWER o el código usado para crear `power_semanal_valle_<año>.tif`? En
particular, se necesitan fecha de descarga, `time-standard`, regla de semana y
método de remuestreo.

## H-015 — Experimentos adicionales de temperatura

**Fecha:** 2026-09-04 — **Tipo:** experimento controlado

Sobre las semanas 2–52 de 2024 se compararon explicaciones alternativas:

| Variante | RMSE | Resultado |
|---|---:|---|
| LST, semana desde 1 de enero | **0.2446 °C** | Mejor receta legítima |
| UTC | 0.2459 °C | No explica la diferencia |
| Semana desde 2 de enero | 0.2845 °C | Empeora |
| Semana desde 3 de enero | 0.3463 °C | Empeora más |
| Desplazamiento optimizado de grilla | 0.2320 °C | Mejora pequeña, sin metadatos |
| Calibración lineal diagnóstica | 0.2372 °C | No resuelve el patrón espacial |

El mejor desplazamiento matemático fue +0.0125° en longitud y +0.025° en
latitud, pero no se adopta: fue escogido contra los mismos datos usados para
evaluarlo, la ganancia es pequeña y no existe evidencia de que esa fuera la
georreferenciación original. Se confirma LST, bloques desde el 1 de enero y
bilineal como reconstrucción transparente.

## D-004 — Reconstrucción provisional completa 2010–2025

**Fecha:** 2026-09-04 — **Estado:** ejecutada, sujeta a información del profesor

Se reconstruyeron los 16 años con T2M bilineal y radiación por vecino cercano.
Se preservaron 32 NetCDF diarios y se generaron 16 GeoTIFF de 104 bandas. Cada
banda tiene cobertura válida en exactamente las 686 celdas aprobadas. Ningún
archivo entregado fue reemplazado.

En temperatura, para S02–S52, la correlación anual varía entre 0.9832 y 0.9959;
el RMSE mediano es 0.3162 °C y el máximo anual 0.4910 °C. Es una concordancia
temporal alta y consistente con diferencias de versión/alineación espacial.


## H-016 — Radiación 2024 no representa toda la historia

**Fecha:** 2026-09-04 — **Tipo:** hallazgo crítico

La réplica exacta por vecino cercano en S02–S52 ocurre en 2018 y 2024, pero no
en la mayoría de años. Entre años, el RMSE de radiación varía entre 0 y 1.0242
MJ/m²/día y la correlación entre 0.6638 y 1. Esto descarta afirmar que el archivo
entregado completo pueda reproducirse con una única descarga actual. Es
compatible con revisiones históricas del producto CERES o con insumos creados
en versiones/fechas distintas, pero la causa final requiere el script o los
crudos del profesor.

**Uso recomendado:** los mapas NASA actuales constituyen una serie completa y
reproducible para modelar; los valores entregados se mantienen como evidencia
de auditoría. Antes del modelo final se comparará el desempeño con y sin
radiación y, como sensibilidad, con radiación bilineal. No se presentará la
serie reconstruida como una réplica exacta del archivo suministrado.

**Archivos:** `04_Experimentos_Temperatura_2024.R`,
`05_Reconstruccion_POWER_2010_2025.R` y
`power_reconstruido_2010_2025/`.

## D-005 — Conservar vecino y bilineal para ambas variables

**Fecha:** 2026-09-06 — **Estado:** aprobada y ejecutada

Por decisión conjunta se generaron las dos versiones de temperatura y las dos
versiones de radiación para 2010–2025. La comparación se hizo únicamente contra
las celdas válidas entregadas por el profesor, nunca contra valores reconstruidos.

Para temperatura en S02–S52, bilineal tuvo menor RMSE en los 16 años: RMSE
mediano anual 0.3162 °C frente a 0.5482 °C de vecino; correlación mediana 0.9930
frente a 0.9794. La evidencia favorece inequívocamente bilineal como candidato
principal, aunque se conserva vecino para sensibilidad.

Para radiación no existe un ganador fuerte: vecino tuvo menor RMSE en 10 años y
bilineal en 6. Las medianas anuales son prácticamente iguales (0.7094 vecino y
0.7098 bilineal), al igual que las correlaciones (0.8438 y 0.8411). Por ello no
se escogerá radiación antes del modelado: ambas entrarán como escenarios y se
compararán también contra un modelo sin radiación.

**Archivos:** `06_Comparacion_Vecino_Bilineal_POWER.R`,
`power_reconstruido_2010_2025/alternativas_vecino_bilineal/` y
`power_reconstruido_2010_2025/comparacion_metodos/`.

## D-006 — Dataset maestro para iniciar el EDA

**Fecha:** 2026-09-06 — **Estado:** construido y verificado

Se creó una tabla larga de 570.752 filas: 686 celdas × 16 años × 52 semanas.
Cada fila identifica celda, coordenadas, altitud, año y semana, e incluye
precipitación CHIRPS y las versiones vecino/bilineal de temperatura y radiación.
Todas las columnas tienen cobertura completa sobre el soporte reconstruido.

El formato principal es RDS para preservar tipos, reducir tamaño y cargarlo
rápidamente. También se guardaron diccionario, control de cobertura y una muestra
CSV de 2.000 filas. Los productos POWER siguen identificados como reconstruidos;
su ausencia de NA no los convierte en mediciones originales a resolución 0.05°.

**Archivo:** `dataset_eda/dataset_eda_2010_2025.rds`.
**Código:** `07_Construir_Dataset_EDA.R`.

## H-017 — Primer EDA sobre el soporte común

**Fecha:** 2026-09-06 — **Tipo:** hallazgo exploratorio

El EDA usa 570.752 filas completas. La precipitación semanal tiene media 50.78
mm, mediana 36.91 mm, desviación 51.86 mm, máximo 660.31 mm y 9.88% de ceros;
por tanto, es asimétrica y con extremos, aspecto que deberá considerarse al
definir la variable/modelo de media.

La resolución efectiva queda cuantificada: por capa, vecino produce normalmente
13–14 valores espaciales distintos de temperatura y 5–6 de radiación; bilineal
produce cerca de 686. Estos últimos son valores interpolados, no observaciones
independientes nuevas.

Las correlaciones cambian con la escala. En datos agrupados espacio-tiempo,
precipitación se correlaciona 0.453 con temperatura bilineal, -0.222 con
radiación bilineal y -0.394 con altitud. En los promedios espaciales son 0.879,
-0.643 y -0.671. La diferencia demuestra que no deben interpretarse relaciones
pooled como efectos causales: mezclan ciclo temporal y gradientes espaciales.

Todas las superficies promedio muestran autocorrelación espacial positiva alta.
Esto es descriptivo, no prueba todavía que el modelo necesite covarianza
espacial: la decisión correcta depende del Moran y semivariograma de los
residuales después de definir la media.

**Archivos:** `08_EDA_Soporte_Comun.R` y `eda_resultados/`.

---

# 10. Plantilla de futuras entradas

```text
Código y fecha:
Tipo: hecho / decisión / propuesta / riesgo / acción
Pregunta:
Alternativas:
Evidencia:
Decisión y justificación:
Consecuencias:
Archivos relacionados:
Estado: pendiente / aprobada / rechazada / sustituida
```

---

# 11. Cumplimiento del enunciado

- Desarrollo en R y metodología principal vista en clase.
- Librerías adicionales solo para gráficos, mapas o carga autorizada.
- Entrega final en un único archivo R autocontenido.
- Informe final máximo 10 páginas, no elaborado en R Markdown.
- Todo resultado acompañado de interpretación contextual.
- Por el uso de asistencia generativa deberá prepararse el documento separado de
  prompts exigido. Esta bitácora no lo sustituye; debe conservarse/exportarse la
  conversación original.
