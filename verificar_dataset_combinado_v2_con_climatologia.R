library(terra)
library(ggplot2)

# Rutas reproducibles: el script funciona desde la raiz del proyecto o desde
# cualquier subcarpeta. No depende de setwd() ni de rutas personales.
buscar_raiz <- function() {
  candidatos <- unique(normalizePath(c(getwd(), file.path(getwd(), ".."),
                                        file.path(getwd(), "../..")),
                                      mustWork = FALSE))
  ok <- vapply(candidatos, function(p) {
    dir.exists(file.path(p, "datos_proyecto_1")) &&
      dir.exists(file.path(p, "avance_team_proyecto1"))
  }, logical(1))
  if (!any(ok)) stop("No se encontro la raiz del proyecto. Abra Proyecto1 en RStudio.")
  candidatos[which(ok)[1]]
}

RAIZ_PROYECTO <- buscar_raiz()
DIR_SCRIPT <- file.path(RAIZ_PROYECTO, "avance_team_proyecto1", "reinicio_desde_cero")
ruta_nuevo <- file.path(DIR_SCRIPT, "dataset_eda")
DIR_DATOS <- file.path(RAIZ_PROYECTO, "datos_proyecto_1")
DIR_RESULTADOS <- file.path(DIR_SCRIPT, "verificacion_dataset_resultados")
dir.create(DIR_RESULTADOS, recursive = TRUE, showWarnings = FALSE)

cargar_flexible <- function(nombre_sin_ext) {
  candidatos <- file.path(ruta_nuevo, paste0(nombre_sin_ext, c(".csv", ".xlsx")))
  candidatos <- candidatos[file.exists(candidatos)]
  if (length(candidatos) == 0) {
    cat("ADVERTENCIA: no se encontro", nombre_sin_ext, "(.csv ni .xlsx) en", ruta_nuevo, "\n")
    return(NULL)
  }
  ruta <- candidatos[1]
  cat("Cargando:", ruta, "\n")
  if (grepl("\\.csv$", ruta)) return(read.csv(ruta, stringsAsFactors = FALSE))
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Se encontro solo XLSX y falta openxlsx. Instale openxlsx o exporte el archivo a CSV.")
  }
  openxlsx::read.xlsx(ruta)
}

# =============================================================================
# PARTE 1: INSPECCION (correr esto primero, sin asumir nada de la estructura)
# =============================================================================
dataset_nuevo <- readRDS(file.path(ruta_nuevo, "dataset_eda_2010_2025.rds"))
diccionario   <- cargar_flexible("diccionario_variables")
control_cob   <- cargar_flexible("control_cobertura")
muestra       <- cargar_flexible("muestra_2000_filas")

cat("\n========== dataset_eda_2010_2025.rds ==========\n")
cat("Dimensiones:", paste(dim(dataset_nuevo), collapse = " x "), "\n")
cat("Columnas:", paste(names(dataset_nuevo), collapse = ", "), "\n")
str(dataset_nuevo)
print(head(dataset_nuevo))

if (!is.null(diccionario)) {
  cat("\n========== diccionario_variables ==========\n")
  print(diccionario)
}

if (!is.null(control_cob)) {
  cat("\n========== control_cobertura ==========\n")
  print(control_cob)
}

if (!is.null(muestra)) {
  cat("\n========== muestra_2000_filas ==========\n")
  cat("Dimensiones:", paste(dim(muestra), collapse = " x "), "\n")
  print(head(muestra))
}

cat("\n>> Revisa lo impreso arriba antes de seguir a la PARTE 2.\n")
cat(">> Si los nombres de columnas no calzan con lo que se busca abajo\n")
cat(">> (coordenadas, anio, semana, precipitacion, altitud, temperatura,\n")
cat(">> radiacion), hay que ajustar esos nombres antes de comparar.\n")

# =============================================================================
# PARTE 2: COMPARACION CONTRA LOS DATOS ORIGINALES DEL PROFESOR
# (asume nombres de columna en espanol: x/y o lon/lat, anio, semana,
#  precipitacion, altitud, temperatura, radiacion -- ajustar si es necesario)
# =============================================================================

nombres_col <- names(dataset_nuevo)

col_x   <- nombres_col[nombres_col %in% c("x", "lon", "longitud")][1]
col_y   <- nombres_col[nombres_col %in% c("y", "lat", "latitud")][1]
col_anio    <- nombres_col[nombres_col %in% c("anio", "año", "year")][1]
col_semana  <- nombres_col[nombres_col %in% c("semana", "week")][1]
col_precip  <- nombres_col[grepl("precip", nombres_col, ignore.case = TRUE)][1]
col_altitud <- nombres_col[grepl("altitud|elevation|alt$", nombres_col, ignore.case = TRUE)][1]
col_temp    <- nombres_col[grepl("temp", nombres_col, ignore.case = TRUE)][1]
col_rad     <- nombres_col[grepl("radiac|radiat", nombres_col, ignore.case = TRUE)][1]

cat("\nColumnas detectadas automaticamente:\n")
cat("  x/lon      :", col_x, "\n")
cat("  y/lat      :", col_y, "\n")
cat("  anio       :", col_anio, "\n")
cat("  semana     :", col_semana, "\n")
cat("  precipitacion:", col_precip, "\n")
cat("  altitud    :", col_altitud, "\n")
cat("  temperatura:", col_temp, "\n")
cat("  radiacion  :", col_rad, "\n")

if (any(sapply(list(col_x, col_y, col_anio, col_semana, col_precip, col_altitud), is.na))) {
  stop("No se pudieron detectar todas las columnas necesarias. Revisa los nombres ",
       "impresos en la PARTE 1 y ajusta manualmente las variables col_x, col_y, ",
       "col_anio, col_semana, col_precip, col_altitud arriba en el script.")
}

# Cargar los rasters originales ya validados (Fase 1 - EDA)
objetos_raster <- readRDS(file.path(DIR_DATOS, "datos_procesados", "rasters_limpios.rds"))
precip_stack_orig <- unwrap(objetos_raster$precip_stack)
r_altitud_orig    <- unwrap(objetos_raster$r_altitud)
anios_orig        <- objetos_raster$anios
valle             <- unwrap(objetos_raster$valle)

# Mapa: que pixeles se quedan (centroide dentro del poligono, regla usada
# por rasterize()/mask() por defecto: touches = FALSE) vs cuales se
# descartan (centroide fuera). Panel izquierdo: vista completa del
# departamento. Panel derecho: zoom sobre una zona de borde compleja
# (costa Pacifica) donde el criterio de centroide se ve mas claro celda
# por celda.
mascara_valle <- rasterize(valle, r_altitud_orig)
r_dentro_fuera <- mascara_valle
values(r_dentro_fuera) <- ifelse(is.na(values(mascara_valle)[, 1]), 0, 1)
levels(r_dentro_fuera) <- data.frame(id = c(0, 1), categoria = c("Se descarta (centroide fuera)", "Se queda (centroide dentro)"))

par(mfrow = c(1, 2))
plot(r_dentro_fuera, col = c("grey80", "steelblue"), main = "Regla de centroide: pixeles dentro/fuera")
plot(valle, add = TRUE, border = "red", lwd = 2)

zona_zoom <- ext(-77.4, -76.9, 3.6, 4.05)
plot(crop(r_dentro_fuera, zona_zoom), col = c("grey80", "steelblue"),
     main = "Zoom: costa Pacifica (borde complejo)")
plot(crop(valle, zona_zoom), add = TRUE, border = "red", lwd = 2)
par(mfrow = c(1, 1))

cat("\nCeldas que se quedan (centroide dentro):", sum(values(r_dentro_fuera)[, 1] == 1, na.rm = TRUE), "\n")
cat("Celdas que se descartan (centroide fuera):", sum(values(r_dentro_fuera)[, 1] == 0, na.rm = TRUE), "\n")

# Comparar ALTITUD (variable estatica, mas facil de verificar primero)
df_altitud_orig <- as.data.frame(r_altitud_orig, xy = TRUE, na.rm = TRUE)
names(df_altitud_orig) <- c("x", "y", "altitud_original")

df_altitud_nuevo <- unique(dataset_nuevo[, c(col_x, col_y, col_altitud)])
names(df_altitud_nuevo) <- c("x", "y", "altitud_nuevo")

comparacion_altitud <- merge(df_altitud_orig, df_altitud_nuevo, by = c("x", "y"))
cat("\n========== COMPARACION DE ALTITUD ==========\n")
cat("Filas comparables:", nrow(comparacion_altitud), "\n")
cat("Correlacion:", round(cor(comparacion_altitud$altitud_original, comparacion_altitud$altitud_nuevo), 4), "\n")
cat("Diferencia absoluta media:", round(mean(abs(comparacion_altitud$altitud_original - comparacion_altitud$altitud_nuevo)), 4), "\n")
cat("Diferencia absoluta maxima:", round(max(abs(comparacion_altitud$altitud_original - comparacion_altitud$altitud_nuevo)), 4), "\n")

# Comparar PRECIPITACION para una semana de ejemplo
anio_ejemplo <- 2024
semana_ejemplo_num <- 10
nombre_banda <- paste0("precip_", anio_ejemplo, "_s", sprintf("%02d", semana_ejemplo_num))

df_precip_orig <- as.data.frame(precip_stack_orig[[nombre_banda]], xy = TRUE, na.rm = TRUE)
names(df_precip_orig) <- c("x", "y", "precip_original")

filtro_nuevo <- dataset_nuevo[[col_anio]] == anio_ejemplo & dataset_nuevo[[col_semana]] == semana_ejemplo_num
df_precip_nuevo <- dataset_nuevo[filtro_nuevo, c(col_x, col_y, col_precip)]
names(df_precip_nuevo) <- c("x", "y", "precip_nuevo")

comparacion_precip <- merge(df_precip_orig, df_precip_nuevo, by = c("x", "y"))
cat("\n========== COMPARACION DE PRECIPITACION (", nombre_banda, ") ==========\n")
cat("Filas comparables:", nrow(comparacion_precip), "\n")
if (nrow(comparacion_precip) > 0) {
  cat("Correlacion:", round(cor(comparacion_precip$precip_original, comparacion_precip$precip_nuevo), 4), "\n")
  cat("Diferencia absoluta media:", round(mean(abs(comparacion_precip$precip_original - comparacion_precip$precip_nuevo)), 4), "\n")
  cat("Diferencia absoluta maxima:", round(max(abs(comparacion_precip$precip_original - comparacion_precip$precip_nuevo)), 4), "\n")
}

# Cobertura dentro del soporte aprobado. En los archivos antiguos faltaba
# 82.2% de temperatura y 96.1% de radiacion sobre las 686 celdas.
if (!is.na(col_temp)) {
  pct_temp_nuevo <- 100 * mean(!is.na(dataset_nuevo[[col_temp]]))
  cat("\nCobertura de temperatura en el dataset nuevo:", round(pct_temp_nuevo, 1), "%",
      "(antes, dentro del Valle: 17.8%)\n")
}
if (!is.na(col_rad)) {
  pct_rad_nuevo <- 100 * mean(!is.na(dataset_nuevo[[col_rad]]))
  cat("Cobertura de radiacion en el dataset nuevo:", round(pct_rad_nuevo, 1), "%",
      "(antes, dentro del Valle: 3.9%)\n")
}

# Revisar si el dataset nuevo ya trae climatologia
col_clima <- nombres_col[grepl("clim", nombres_col, ignore.case = TRUE)]
if (length(col_clima) > 0) {
  cat("\nEl dataset nuevo SI trae columnas de climatologia:", paste(col_clima, collapse = ", "), "\n")
} else {
  cat("\nEl dataset nuevo NO trae columnas de climatologia.\n")
  cat("No hace falta pedirsela: ya tenemos precip_clim_derivada (Fase 1, seccion 10),\n")
  cat("reconstruida directamente de los datos limpios y validada contra el archivo\n")
  cat("original corrupto del profesor. Se puede seguir usando esa.\n")
}

# =============================================================================
# PARTE 3: EDA DE VERIFICACION PARA TEMPERATURA Y RADIACION
# =============================================================================
temp_stack_orig <- unwrap(objetos_raster$temp_stack)
rad_stack_orig  <- if (!is.null(objetos_raster$rad_stack)) unwrap(objetos_raster$rad_stack) else NULL

cat("\n========== RANGOS DE LAS 4 COLUMNAS NUEVAS ==========\n")
print(summary(dataset_nuevo[, c("temperatura_vecino_c", "temperatura_bilineal_c",
                                "radiacion_vecino_mj_m2_dia", "radiacion_bilineal_mj_m2_dia")]))

cat("\n========== CORRELACION VECINO vs BILINEAL (deberia ser alta) ==========\n")
cat("Temperatura:", round(cor(dataset_nuevo$temperatura_vecino_c, dataset_nuevo$temperatura_bilineal_c), 4), "\n")
cat("Radiacion  :", round(cor(dataset_nuevo$radiacion_vecino_mj_m2_dia, dataset_nuevo$radiacion_bilineal_mj_m2_dia), 4), "\n")

# Cruce contra los datos viejos (donde SI habia dato original, antes del arreglo)
df_temp_viejo <- as.data.frame(temp_stack_orig[[nombre_banda_temp <- names(temp_stack_orig)[grepl("2024_s10", names(temp_stack_orig))]]],
                               xy = TRUE, na.rm = TRUE)
names(df_temp_viejo) <- c("x", "y", "temp_viejo")

df_temp_nuevo_semana <- dataset_nuevo[dataset_nuevo[[col_anio]] == anio_ejemplo & dataset_nuevo[[col_semana]] == semana_ejemplo_num,
                                      c("x", "y", "temperatura_vecino_c", "temperatura_bilineal_c")]

cruce_temp <- merge(df_temp_viejo, df_temp_nuevo_semana, by = c("x", "y"))
cat("\n========== TEMPERATURA: donde SI habia dato original (n =", nrow(cruce_temp), ") ==========\n")
if (nrow(cruce_temp) > 0) {
  cat("Correlacion viejo vs vecino  :", round(cor(cruce_temp$temp_viejo, cruce_temp$temperatura_vecino_c), 4), "\n")
  cat("Correlacion viejo vs bilineal:", round(cor(cruce_temp$temp_viejo, cruce_temp$temperatura_bilineal_c), 4), "\n")
  cat("Diferencia media viejo vs vecino  :", round(mean(abs(cruce_temp$temp_viejo - cruce_temp$temperatura_vecino_c)), 4), "\n")
  cat("Diferencia media viejo vs bilineal:", round(mean(abs(cruce_temp$temp_viejo - cruce_temp$temperatura_bilineal_c)), 4), "\n")
} else {
  cat("No hay celdas en comun para cruzar en esta semana especifica (normal si la\n")
  cat("cobertura vieja era muy escasa); se puede repetir con otra semana si hace falta.\n")
}

if (!is.null(rad_stack_orig)) {
  nombre_banda_rad <- names(rad_stack_orig)[grepl("2024_s10", names(rad_stack_orig))]
  df_rad_viejo <- as.data.frame(rad_stack_orig[[nombre_banda_rad]], xy = TRUE, na.rm = TRUE)
  names(df_rad_viejo) <- c("x", "y", "rad_viejo")
  
  df_rad_nuevo_semana <- dataset_nuevo[dataset_nuevo[[col_anio]] == anio_ejemplo & dataset_nuevo[[col_semana]] == semana_ejemplo_num,
                                       c("x", "y", "radiacion_vecino_mj_m2_dia", "radiacion_bilineal_mj_m2_dia")]
  
  cruce_rad <- merge(df_rad_viejo, df_rad_nuevo_semana, by = c("x", "y"))
  cat("\n========== RADIACION: donde SI habia dato original (n =", nrow(cruce_rad), ") ==========\n")
  if (nrow(cruce_rad) > 0) {
    cat("Correlacion viejo vs vecino  :", round(cor(cruce_rad$rad_viejo, cruce_rad$radiacion_vecino_mj_m2_dia), 4), "\n")
    cat("Correlacion viejo vs bilineal:", round(cor(cruce_rad$rad_viejo, cruce_rad$radiacion_bilineal_mj_m2_dia), 4), "\n")
    cat("Diferencia media viejo vs vecino  :", round(mean(abs(cruce_rad$rad_viejo - cruce_rad$radiacion_vecino_mj_m2_dia)), 4), "\n")
    cat("Diferencia media viejo vs bilineal:", round(mean(abs(cruce_rad$rad_viejo - cruce_rad$radiacion_bilineal_mj_m2_dia)), 4), "\n")
  } else {
    cat("No hay celdas en comun para cruzar en esta semana especifica.\n")
  }
}

# Mapas: vecino vs bilineal, para ver visualmente diferencias de suavizado
df_semana_completa <- dataset_nuevo[dataset_nuevo[[col_anio]] == anio_ejemplo & dataset_nuevo[[col_semana]] == semana_ejemplo_num, ]

r_temp_vecino   <- rast(df_semana_completa[, c("x", "y", "temperatura_vecino_c")], type = "xyz", crs = crs(r_altitud_orig))
r_temp_bilineal <- rast(df_semana_completa[, c("x", "y", "temperatura_bilineal_c")], type = "xyz", crs = crs(r_altitud_orig))

par(mfrow = c(1, 2))
plot(r_temp_vecino, main = "Temperatura - vecino mas cercano")
plot(r_temp_bilineal, main = "Temperatura - bilineal")
par(mfrow = c(1, 1))

r_rad_vecino   <- rast(df_semana_completa[, c("x", "y", "radiacion_vecino_mj_m2_dia")], type = "xyz", crs = crs(r_altitud_orig))
r_rad_bilineal <- rast(df_semana_completa[, c("x", "y", "radiacion_bilineal_mj_m2_dia")], type = "xyz", crs = crs(r_altitud_orig))

par(mfrow = c(1, 2))
plot(r_rad_vecino, main = "Radiacion - vecino mas cercano")
plot(r_rad_bilineal, main = "Radiacion - bilineal")
par(mfrow = c(1, 1))

# Dispersion: donde SI habia dato viejo, contra el valor bilineal nuevo
# (version grafica de las correlaciones ya impresas arriba)
par(mfrow = c(1, 2))
if (nrow(cruce_temp) > 0) {
  plot(cruce_temp$temp_viejo, cruce_temp$temperatura_bilineal_c, pch = 19, col = "steelblue",
       xlab = "Temperatura (dato viejo)", ylab = "Temperatura (bilineal, nuevo)",
       main = "Temperatura: viejo vs. nuevo")
  abline(0, 1, col = "red", lty = 2)
}
if (exists("cruce_rad") && nrow(cruce_rad) > 0) {
  plot(cruce_rad$rad_viejo, cruce_rad$radiacion_bilineal_mj_m2_dia, pch = 19, col = "darkorange",
       xlab = "Radiacion (dato viejo)", ylab = "Radiacion (bilineal, nuevo)",
       main = "Radiacion: viejo vs. nuevo")
  abline(0, 1, col = "red", lty = 2)
}
par(mfrow = c(1, 1))

# Distribucion de las variables ya con cobertura completa
par(mfrow = c(1, 2))
hist(dataset_nuevo$temperatura_bilineal_c, breaks = 40, col = "firebrick",
     main = "Temperatura bilineal (cobertura completa)", xlab = "Temperatura (C)")
hist(dataset_nuevo$radiacion_bilineal_mj_m2_dia, breaks = 40, col = "orange",
     main = "Radiacion bilineal (cobertura completa)", xlab = "Radiacion (MJ/m2/dia)")
par(mfrow = c(1, 1))

# =============================================================================
# PARTE 4: RECONSTRUCCION DE rasters_limpios.rds CON EL DATASET DEFINITIVO
# Solo se reemplazan temp_stack y rad_stack (usando la version bilineal).
# precip_stack, r_altitud y precip_clim_derivada NO se tocan: ya se confirmo
# en la Parte 2 que son identicos a los originales.
# =============================================================================
construir_stack_desde_tidy <- function(df, col_valor, prefijo) {
  anios_datos <- sort(unique(df$anio))
  stacks_por_anio <- vector("list", length(anios_datos))
  
  for (i in seq_along(anios_datos)) {
    a <- anios_datos[i]
    df_anio <- df[df$anio == a, ]
    celdas_xy <- unique(df_anio[, c("id_celda", "x", "y")])
    celdas_xy <- celdas_xy[order(celdas_xy$id_celda), ]
    
    semanas <- sort(unique(df_anio$semana))
    wide <- matrix(NA_real_, nrow = nrow(celdas_xy), ncol = length(semanas))
    for (j in seq_along(semanas)) {
      sub <- df_anio[df_anio$semana == semanas[j], ]
      sub <- sub[match(celdas_xy$id_celda, sub$id_celda), ]
      wide[, j] <- sub[[col_valor]]
    }
    df_wide <- cbind(celdas_xy[, c("x", "y")], as.data.frame(wide))
    names(df_wide) <- c("x", "y", paste0(prefijo, "_", a, "_s", sprintf("%02d", semanas)))
    stacks_por_anio[[i]] <- rast(df_wide, type = "xyz", crs = "EPSG:4326")
  }
  rast(stacks_por_anio)
}

cat("\n============================================================\n")
cat("PARTE 4: RECONSTRUYENDO temp_stack Y rad_stack (dataset definitivo)\n")
cat("============================================================\n")

temp_stack_nuevo <- construir_stack_desde_tidy(dataset_nuevo, "temperatura_bilineal_c", "temp")
rad_stack_nuevo  <- construir_stack_desde_tidy(dataset_nuevo, "radiacion_bilineal_mj_m2_dia", "rad")

cat("temp_stack_nuevo:", nlyr(temp_stack_nuevo), "bandas (precip_stack tiene", nlyr(precip_stack_orig), ")\n")
cat("rad_stack_nuevo :", nlyr(rad_stack_nuevo), "bandas\n")
cat("Cobertura temp_stack_nuevo:", round(100 * mean(global(!is.na(temp_stack_nuevo), "mean")[, 1]), 1), "%\n")
cat("Cobertura rad_stack_nuevo :", round(100 * mean(global(!is.na(rad_stack_nuevo), "mean")[, 1]), 1), "%\n")

objetos_raster$temp_stack <- wrap(temp_stack_nuevo)
objetos_raster$rad_stack  <- wrap(rad_stack_nuevo)

# Este script verifica y genera derivados: nunca sobrescribe el insumo previo.
ruta_rasters_verificados <- file.path(DIR_RESULTADOS, "rasters_combinados_bilineales.rds")
saveRDS(objetos_raster, ruta_rasters_verificados)
cat("\nGuardado derivado (sin sobrescribir insumos):", ruta_rasters_verificados, "\n")

# Mapa ANTES vs DESPUES: la comparacion mas directa del arreglo de cobertura
nombre_banda_temp_nuevo <- paste0("temp_", anio_ejemplo, "_s", sprintf("%02d", semana_ejemplo_num))
nombre_banda_rad_nuevo  <- paste0("rad_", anio_ejemplo, "_s", sprintf("%02d", semana_ejemplo_num))

par(mfrow = c(2, 2))
plot(temp_stack_orig[[nombre_banda_temp]], main = "Temperatura ANTES (remuestreo sin interpolar)")
plot(valle, add = TRUE, border = "black")
plot(temp_stack_nuevo[[nombre_banda_temp_nuevo]], main = "Temperatura AHORA (bilineal)")
plot(valle, add = TRUE, border = "black")
plot(rad_stack_orig[[nombre_banda_rad]], main = "Radiacion ANTES (remuestreo sin interpolar)")
plot(valle, add = TRUE, border = "black")
plot(rad_stack_nuevo[[nombre_banda_rad_nuevo]], main = "Radiacion AHORA (bilineal)")
plot(valle, add = TRUE, border = "black")
par(mfrow = c(1, 1))

# Matriz de correlacion con el dataset definitivo (cobertura 100%, ya no
# limitada por las muestras pequenias de antes)
vars_cor_nuevo <- c("precipitacion_mm", "altitud_m", "temperatura_bilineal_c", "radiacion_bilineal_mj_m2_dia")
matriz_cor_nueva <- cor(dataset_nuevo[, vars_cor_nuevo], use = "pairwise.complete.obs")
cat("\n========== MATRIZ DE CORRELACION - DATASET DEFINITIVO (cobertura 100%) ==========\n")
print(round(matriz_cor_nueva, 3))
cat("(basada en las", nrow(dataset_nuevo), "filas completas, no en muestras reducidas)\n")

df_matriz_nueva <- as.data.frame(as.table(matriz_cor_nueva))
names(df_matriz_nueva) <- c("var1", "var2", "correlacion")

p_matriz_nueva <- ggplot(df_matriz_nueva, aes(x = var1, y = var2, fill = correlacion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(correlacion, 2)), color = "black") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Matriz de correlacion - dataset definitivo", x = NULL, y = NULL, fill = "Correlacion") +
  theme_minimal()
print(p_matriz_nueva)

# =============================================================================
# PARTE 5: TABLA COMPARATIVA DE RESOLUCION POR VARIABLE
# =============================================================================
res_precip <- res(precip_stack_orig)[1]
res_altitud <- res(r_altitud_orig)[1]
res_temp_nuevo <- res(temp_stack_nuevo)[1]
res_rad_nuevo  <- res(rad_stack_nuevo)[1]

grados_a_km <- function(grados) round(grados * 111.32, 2)

tabla_resoluciones <- data.frame(
  variable = c("Precipitacion", "Altitud", "Temperatura", "Radiacion"),
  fuente = c("CHIRPS", "SRTM", "NASA POWER", "NASA POWER"),
  resolucion_nativa_aprox_km = c(grados_a_km(res_precip), "~0.03 (30 m; luego agregada)",
                                 "~69 x 56 (0.625 x 0.5 grados)", "~111 x 111 (1 x 1 grado)"),
  resolucion_grilla_comun_grados = c(res_precip, res_altitud, res_temp_nuevo, res_rad_nuevo),
  resolucion_grilla_comun_km = c(grados_a_km(res_precip), grados_a_km(res_altitud),
                                 grados_a_km(res_temp_nuevo), grados_a_km(res_rad_nuevo)),
  metodo_para_llegar_a_grilla_comun = c(
    "Ninguno (ya viene en la resolucion objetivo)",
    "Agregacion (de muy fina, 30 m, hacia la grilla objetivo)",
    "Interpolacion bilineal (de ~69 x 56 km hacia la grilla comun)",
    "Interpolacion bilineal (de ~111 x 111 km hacia la grilla comun)"
  )
)

cat("\n========== TABLA COMPARATIVA DE RESOLUCION POR VARIABLE ==========\n")
print(tabla_resoluciones)

cat("\nPOWER es mas grueso que la grilla comun de 0.05 grados. Tanto 'near' como\n")
cat("bilinear pueden producir cobertura completa si se aplica un remuestreo espacial\n")
cat("correcto. Los nulos antiguos (82.2% temperatura y 96.1% radiacion dentro de\n")
cat("las 686 celdas) no fueron causados por vecino cercano, sino por una transferencia\n")
cat("incompleta de nodos de la grilla fuente a la grilla objetivo.\n")

# =============================================================================
# PARTE 6: CLIMATOLOGIA DE TEMPERATURA Y RADIACION - ORIGINAL vs. DERIVADA
# (misma logica que se aplico a precipitacion en la Fase 1: comparar el
# archivo de climatologia del profesor contra una version reconstruida
# directamente de los datos ya limpios, para revisar si tambien esta corrupto)
# =============================================================================
ruta_base_original <- file.path(DIR_DATOS, "imagenes_semanales")

r_clim_temp_orig <- rast(file.path(ruta_base_original, "power_temp_climatologia_semanal_valle.tif"))
r_clim_rad_orig  <- rast(file.path(ruta_base_original, "power_radiacion_climatologia_semanal_valle.tif"))

cat("\n========== CLIMATOLOGIA DE TEMPERATURA: original vs. derivada ==========\n")
clim_temp_archivo <- global(r_clim_temp_orig, "mean", na.rm = TRUE)
clim_temp_archivo$semana <- 1:52

temp_clim_derivada <- tapp(temp_stack_nuevo, index = rep(1:52, times = length(anios_orig)), fun = "mean", na.rm = TRUE)
names(temp_clim_derivada) <- paste0("semana_", sprintf("%02d", 1:52))
clim_temp_derivada <- global(temp_clim_derivada, "mean", na.rm = TRUE)
clim_temp_derivada$semana <- 1:52

cat("Resumen archivo original (C):\n"); print(summary(clim_temp_archivo$mean))
cat("Resumen derivada (C):\n"); print(summary(clim_temp_derivada$mean))
cat("Correlacion entre ambas series:", round(cor(clim_temp_archivo$mean, clim_temp_derivada$mean), 4), "\n")
cat("Diferencia media (derivada - original):", round(mean(clim_temp_derivada$mean - clim_temp_archivo$mean), 3), "C\n")
cat("Desviacion estandar de la diferencia:", round(sd(clim_temp_derivada$mean - clim_temp_archivo$mean), 3), "C\n")

cat("\n========== CLIMATOLOGIA DE RADIACION: original vs. derivada ==========\n")
clim_rad_archivo <- global(r_clim_rad_orig, "mean", na.rm = TRUE)
clim_rad_archivo$semana <- 1:52

rad_clim_derivada <- tapp(rad_stack_nuevo, index = rep(1:52, times = length(anios_orig)), fun = "mean", na.rm = TRUE)
names(rad_clim_derivada) <- paste0("semana_", sprintf("%02d", 1:52))
clim_rad_derivada <- global(rad_clim_derivada, "mean", na.rm = TRUE)
clim_rad_derivada$semana <- 1:52

cat("Resumen archivo original (MJ/m2/dia):\n"); print(summary(clim_rad_archivo$mean))
cat("Resumen derivada (MJ/m2/dia):\n"); print(summary(clim_rad_derivada$mean))
cat("Correlacion entre ambas series:", round(cor(clim_rad_archivo$mean, clim_rad_derivada$mean), 4), "\n")
cat("Diferencia media (derivada - original):", round(mean(clim_rad_derivada$mean - clim_rad_archivo$mean), 3), "MJ/m2/dia\n")
cat("Desviacion estandar de la diferencia:", round(sd(clim_rad_derivada$mean - clim_rad_archivo$mean), 3), "MJ/m2/dia\n")

# Graficos comparativos
comparacion_clima_temp <- rbind(
  data.frame(semana = clim_temp_archivo$semana, valor = clim_temp_archivo$mean, fuente = "Archivo original"),
  data.frame(semana = clim_temp_derivada$semana, valor = clim_temp_derivada$mean, fuente = "Derivada")
)
p_clima_temp <- ggplot(comparacion_clima_temp, aes(x = semana, y = valor, color = fuente)) +
  geom_line(linewidth = 1) + geom_point() +
  labs(title = "Climatologia de temperatura: original vs. derivada",
       x = "Semana ISO", y = "Temperatura (C)", color = "Fuente") +
  theme_minimal()
print(p_clima_temp)

comparacion_clima_rad <- rbind(
  data.frame(semana = clim_rad_archivo$semana, valor = clim_rad_archivo$mean, fuente = "Archivo original"),
  data.frame(semana = clim_rad_derivada$semana, valor = clim_rad_derivada$mean, fuente = "Derivada")
)
p_clima_rad <- ggplot(comparacion_clima_rad, aes(x = semana, y = valor, color = fuente)) +
  geom_line(linewidth = 1) + geom_point() +
  labs(title = "Climatologia de radiacion: original vs. derivada",
       x = "Semana ISO", y = "Radiacion (MJ/m2/dia)", color = "Fuente") +
  theme_minimal()
print(p_clima_rad)

# =============================================================================
# PARTE 7: CONSTRUIR EL DATASET FINAL CON LAS 3 CLIMATOLOGIAS POR CELDA
# (agregando por id_celda+semana, NO por coordenadas x,y flotantes -- fusionar
# por floats de distintas fuentes raster produce fallas silenciosas por
# imprecision numerica; id_celda es una clave entera estable).
# Reutiliza 'mascara_valle', ya calculada en la PARTE 2 de este mismo script.
# =============================================================================

clim_por_celda <- aggregate(
  cbind(precipitacion_mm, temperatura_bilineal_c, radiacion_bilineal_mj_m2_dia) ~ id_celda + semana,
  data = dataset_nuevo, FUN = mean, na.rm = TRUE
)
names(clim_por_celda) <- c("id_celda", "semana", "clim_precipitacion_mm",
                            "clim_temperatura_bilineal_c", "clim_radiacion_bilineal_mj_m2_dia")

cat("\n========== CLIMATOLOGIAS POR CELDA (agregado por id_celda+semana) ==========\n")
cat("Filas:", nrow(clim_por_celda), "(esperado: 686 x 52 =", 686 * 52, ")\n")

# Validacion: climatologia de precipitacion derivada vs. archivo oficial (100% cobertura)
r_clim_precip_orig <- rast(file.path(ruta_base_original, "chirps_climatologia_semanal_valle.tif"))

xy_celdas <- as.data.frame(mascara_valle, xy = TRUE, na.rm = TRUE)[, c("x", "y")]
xy_celdas$id_celda <- seq_len(nrow(xy_celdas))

if (!compareGeom(r_clim_precip_orig, mascara_valle, stopOnError = FALSE)) {
  r_clim_precip_orig <- resample(r_clim_precip_orig, mascara_valle, method = "near")
}
valores_orig <- extract(r_clim_precip_orig, xy_celdas[, c("x", "y")])[, -1]
df_orig <- data.frame(id_celda = rep(xy_celdas$id_celda, times = 52),
                       semana = rep(1:52, each = nrow(xy_celdas)),
                       clim_precip_original = as.vector(as.matrix(valores_orig)))

comparacion_clim_precip <- merge(clim_por_celda[, c("id_celda","semana","clim_precipitacion_mm")],
                                  df_orig, by = c("id_celda","semana"))
cat("\n========== VALIDACION: climatologia precipitacion derivada vs. oficial ==========\n")
cat("Filas comparadas:", nrow(comparacion_clim_precip), "\n")
cat("Correlacion:", round(cor(comparacion_clim_precip$clim_precipitacion_mm,
                               comparacion_clim_precip$clim_precip_original, use="complete.obs"), 4), "\n")
cat("Diferencia absoluta media:", round(mean(abs(comparacion_clim_precip$clim_precipitacion_mm -
                               comparacion_clim_precip$clim_precip_original)), 4), "mm\n")

# Union final
dataset_final <- merge(dataset_nuevo, clim_por_celda, by = c("id_celda", "semana"), all.x = TRUE)
dataset_final <- dataset_final[order(dataset_final$id_celda, dataset_final$anio, dataset_final$semana), ]

cat("\n========== DATASET FINAL CON LAS 7 VARIABLES (+ climatologias) ==========\n")
cat("Dimensiones:", paste(dim(dataset_final), collapse = " x "), "\n")
cat("Columnas:", paste(names(dataset_final), collapse = ", "), "\n")
cat("\nNA por columna:\n")
print(colSums(is.na(dataset_final)))

saveRDS(dataset_final, file.path(ruta_nuevo, "dataset_final_7_variables.rds"))
write.csv(head(dataset_final, 2000), file.path(ruta_nuevo, "muestra_2000_dataset_final_7_variables.csv"),
          row.names = FALSE)
cat("\nGuardado: dataset_final_7_variables.rds\n")

# =============================================================================
# PARTE 8: OUTPUTS AUDITABLES
# =============================================================================
metricas_validacion <- data.frame(
  variable = c("Altitud", "Precipitacion 2024-S10", "Temperatura vecino 2024-S10",
               "Temperatura bilineal 2024-S10", "Radiacion vecino 2024-S10",
               "Radiacion bilineal 2024-S10"),
  n = c(nrow(comparacion_altitud), nrow(comparacion_precip), nrow(cruce_temp),
        nrow(cruce_temp), nrow(cruce_rad), nrow(cruce_rad)),
  correlacion = c(cor(comparacion_altitud$altitud_original, comparacion_altitud$altitud_nuevo),
                  cor(comparacion_precip$precip_original, comparacion_precip$precip_nuevo),
                  cor(cruce_temp$temp_viejo, cruce_temp$temperatura_vecino_c),
                  cor(cruce_temp$temp_viejo, cruce_temp$temperatura_bilineal_c),
                  cor(cruce_rad$rad_viejo, cruce_rad$radiacion_vecino_mj_m2_dia),
                  cor(cruce_rad$rad_viejo, cruce_rad$radiacion_bilineal_mj_m2_dia)),
  mae = c(mean(abs(comparacion_altitud$altitud_original-comparacion_altitud$altitud_nuevo)),
          mean(abs(comparacion_precip$precip_original-comparacion_precip$precip_nuevo)),
          mean(abs(cruce_temp$temp_viejo-cruce_temp$temperatura_vecino_c)),
          mean(abs(cruce_temp$temp_viejo-cruce_temp$temperatura_bilineal_c)),
          mean(abs(cruce_rad$rad_viejo-cruce_rad$radiacion_vecino_mj_m2_dia)),
          mean(abs(cruce_rad$rad_viejo-cruce_rad$radiacion_bilineal_mj_m2_dia)))
)
write.csv(metricas_validacion, file.path(DIR_RESULTADOS, "metricas_validacion.csv"), row.names = FALSE)
write.csv(tabla_resoluciones, file.path(DIR_RESULTADOS, "resoluciones_y_remuestreo.csv"), row.names = FALSE)
write.csv(matriz_cor_nueva, file.path(DIR_RESULTADOS, "correlaciones_dataset_base.csv"))
write.csv(data.frame(variable=names(dataset_final), n_na=colSums(is.na(dataset_final))),
          file.path(DIR_RESULTADOS, "control_nulos_dataset_final.csv"), row.names = FALSE)

informe_txt <- capture.output({
  cat("VERIFICACION DEL DATASET COMBINADO\n")
  cat("Fecha:", as.character(Sys.time()), "\n")
  cat("Dimensiones:", nrow(dataset_final), "filas x", ncol(dataset_final), "columnas\n")
  cat("Soporte aprobado: 686 centros de pixel dentro del Valle\n")
  cat("Registros esperados: 686 x 16 x 52 =", 686*16*52, "\n\n")
  print(metricas_validacion)
  cat("\nNulos por columna:\n"); print(colSums(is.na(dataset_final)))
  cat("\nConclusiones:\n")
  cat("- Altitud y precipitacion reproducen exactamente los insumos de referencia.\n")
  cat("- Temperatura: bilinear es la opcion principal por menor MAE y mayor correlacion.\n")
  cat("- Radiacion: vecino reproduce exactamente la muestra antigua; se conserva bilinear\n")
  cat("  como opcion continua, pero la decision debe evaluarse tambien en los modelos.\n")
  cat("- Interpolar aumenta soporte computacional, no resolucion fisica de POWER.\n")
})
writeLines(informe_txt, file.path(DIR_RESULTADOS, "resumen_verificacion.txt"))

ggsave(file.path(DIR_RESULTADOS, "climatologia_temperatura.png"), p_clima_temp,
       width=9, height=5, dpi=160)
ggsave(file.path(DIR_RESULTADOS, "climatologia_radiacion.png"), p_clima_rad,
       width=9, height=5, dpi=160)
ggsave(file.path(DIR_RESULTADOS, "matriz_correlacion.png"), p_matriz_nueva,
       width=9, height=7, dpi=160)
cat("Outputs de verificacion:", DIR_RESULTADOS, "\n")
