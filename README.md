# Proyecto 1 — entrega limpia hasta EDA descriptivo

Esta carpeta tiene **dos scripts principales para presentar y ejecutar**, en
este orden:

1. `01_Preparar_Verificar_Dataset.R`
   - valida el soporte de 686 centros de pixel;
   - verifica altitud y precipitacion contra los insumos de referencia;
   - compara vecino y bilinear para temperatura y radiacion;
   - calcula las tres climatologias por celda y banda semanal;
   - construye `dataset_eda/dataset_final_7_variables.rds`;
   - guarda evidencia en `verificacion_dataset_resultados/`.
2. `EDA_final.R`
   - realiza exclusivamente el EDA descriptivo no espacial;
   - muestra estadisticas y graficos en RStudio;
   - guarda tablas y figuras en `eda_final_resultados/`.

El punto actual del proyecto termina al ejecutar `EDA_final.R`. El EDA espacial,
Moran, semivariogramas y modelos no forman parte todavía de esta entrega.

La siguiente fase comienza con `02_00_Cierre_Decisiones_Espaciales.R`. Este
archivo solo audita y fija las decisiones de entrada; todavía no calcula Moran,
variogramas ni kriging.

## Cómo ejecutar

Abra la carpeta raíz `Proyecto1` en RStudio. No hace falta usar `setwd()`:

```r
source("avance_team_proyecto1/reinicio_desde_cero/01_Preparar_Verificar_Dataset.R")
source("avance_team_proyecto1/reinicio_desde_cero/EDA_final.R")
```

## Organización

- `_soporte_reproducibilidad/`: módulos que explican cómo se auditó y construyó
  el insumo base. Son necesarios para trazabilidad, pero no son puntos de
  entrada de la presentación.
- `_archivo_historico/`: pilotos, experimentos y versiones reemplazadas. No se
  ejecutan; se conservan para poder recuperar decisiones anteriores.
- `dataset_eda/`: dataset base y dataset final.
- `verificacion_dataset_resultados/`: evidencia de construcción y validación.
- `eda_final_resultados/`: resultados del EDA descriptivo final.
- `BITACORA_PROYECTO.md`: decisiones cronológicas, riesgos y pendientes.

## Variables principales

La respuesta es `precipitacion_mm`. Las covariables son altitud, temperatura,
radiación y sus climatologías semanales. Se conservan las versiones vecino y
bilinear de POWER para trazabilidad, aunque temperatura bilinear es la candidata
principal. La variante de radiación se decidirá mediante validación de modelos.

## Importante

Interpolar POWER a la grilla de 0,05 grados unifica el soporte computacional,
pero no crea observaciones meteorológicas nuevas ni mejora su resolución física
nativa. Los extremos de Tukey se reportan y no se eliminan automáticamente.
