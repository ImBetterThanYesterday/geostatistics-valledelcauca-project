# =============================================================================
# FASE 3 - EDA SOBRE EL SOPORTE ESPACIAL-TEMPORAL COMUN
#
# Entrada: dataset maestro (686 celdas x 16 años x 52 semanas).
# Salidas: consola, panel Plots, PDF multipagina, CSV y objetos resumen.
# No modifica los datos de entrada ni decide todavía el modelo final.
# =============================================================================

library(terra)
library(ggplot2)
library(gridExtra)

detectar_raiz <- function() {
  p <- getwd()
  if (requireNamespace("rstudioapi",quietly=TRUE) && rstudioapi::isAvailable()) {
    a <- tryCatch(rstudioapi::getActiveDocumentContext()$path,error=function(e) "")
    if (nzchar(a)) p <- dirname(a)
  }
  repeat {
    if (dir.exists(file.path(p,"datos_proyecto_1"))) return(normalizePath(p))
    q <- dirname(p); if (identical(p,q)) break; p <- q
  }
  stop("No se encontro la raiz. Abra Proyecto1.Rproj.")
}

RAIZ <- detectar_raiz()
BASE <- file.path(RAIZ,"avance_team_proyecto1","reinicio_desde_cero")
DIR_IN <- file.path(BASE,"dataset_eda")
DIR_OUT <- file.path(BASE,"eda_resultados")
dir.create(DIR_OUT,recursive=TRUE,showWarnings=FALSE)
source(file.path(RAIZ,"avance_team_proyecto1","funciones_auxiliares.R"))

datos <- readRDS(file.path(DIR_IN,"dataset_eda_2010_2025.rds"))
s <- readRDS(file.path(BASE,"soporte_resultados","soporte_maestro.rds"))
mascara <- rast(s$mascara); valle <- vect(s$valle); celdas <- s$celdas

variables <- c("precipitacion_mm","temperatura_vecino_c",
               "temperatura_bilineal_c","radiacion_vecino_mj_m2_dia",
               "radiacion_bilineal_mj_m2_dia","altitud_m")

resumen_variable <- function(x,nombre) {
  q <- quantile(x,c(.25,.5,.75,.95,.99),na.rm=TRUE)
  data.frame(variable=nombre,n=sum(is.finite(x)),na=sum(!is.finite(x)),
    media=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE),min=min(x,na.rm=TRUE),
    p25=q[1],mediana=q[2],p75=q[3],p95=q[4],p99=q[5],max=max(x,na.rm=TRUE),
    pct_ceros=100*mean(x==0,na.rm=TRUE))
}
tabla_resumen <- do.call(rbind,lapply(variables,function(v)
  resumen_variable(datos[[v]],v)))
rownames(tabla_resumen) <- NULL

# Serie regional y climatologia semanal. La media es sobre las mismas 686
# celdas en cada fecha, por lo que no cambia por disponibilidad espacial.
serie <- aggregate(datos[,variables[1:5]],
                   by=list(anio=datos$anio,semana=datos$semana),FUN=mean)
serie$fecha <- as.Date(paste0(serie$anio,"-01-01"))+(serie$semana-1)*7
climatologia <- aggregate(datos[,variables[1:5]],by=list(semana=datos$semana),FUN=mean)
anual <- aggregate(datos[,variables[1:5]],by=list(anio=datos$anio),FUN=mean)

# Promedios espaciales de largo plazo: exactamente una fila por celda.
espacial <- aggregate(datos[,variables],
  by=list(id_celda=datos$id_celda,celda_raster=datos$celda_raster,
          x=datos$x,y=datos$y),FUN=mean)

# Tres correlaciones con interpretaciones diferentes.
cor_pooled <- cor(datos[,variables],use="complete.obs")
cor_espacial <- cor(espacial[,variables],use="complete.obs")

# Anomalias respecto a la climatologia de cada semana: remueven el ciclo anual
# antes de estudiar co-movimiento entre años/celdas.
medias_sem <- aggregate(datos[,variables[1:5]],list(semana=datos$semana),mean)
idx <- match(datos$semana,medias_sem$semana)
anomalias <- datos[,variables[1:5]]-medias_sem[idx,variables[1:5]]
cor_anomalias <- cor(anomalias,use="complete.obs")

matriz_a_largo <- function(m,tipo) {
  z <- as.data.frame(as.table(m)); names(z) <- c("variable_1","variable_2","correlacion")
  z$tipo <- tipo; z
}
correlaciones <- rbind(matriz_a_largo(cor_pooled,"pooled_espacio_tiempo"),
                       matriz_a_largo(cor_espacial,"promedios_espaciales"),
                       matriz_a_largo(cor_anomalias,"anomalias_semanales"))

# Diferencia directa entre métodos sobre exactamente las mismas observaciones.
comparacion_metodos <- data.frame(
  variable=c("temperatura","radiacion"),
  MAE_vecino_vs_bilineal=c(
    mean(abs(datos$temperatura_vecino_c-datos$temperatura_bilineal_c)),
    mean(abs(datos$radiacion_vecino_mj_m2_dia-datos$radiacion_bilineal_mj_m2_dia))),
  RMSE_vecino_vs_bilineal=c(
    sqrt(mean((datos$temperatura_vecino_c-datos$temperatura_bilineal_c)^2)),
    sqrt(mean((datos$radiacion_vecino_mj_m2_dia-datos$radiacion_bilineal_mj_m2_dia)^2))),
  correlacion=c(cor(datos$temperatura_vecino_c,datos$temperatura_bilineal_c),
                cor(datos$radiacion_vecino_mj_m2_dia,
                    datos$radiacion_bilineal_mj_m2_dia)))

# Resolucion efectiva: cantidad mediana de valores distintos por semana-año.
unicos <- do.call(rbind,lapply(split(datos,interaction(datos$anio,datos$semana)),function(z)
  data.frame(anio=z$anio[1],semana=z$semana[1],
    temp_vecino=length(unique(round(z$temperatura_vecino_c,6))),
    temp_bilineal=length(unique(round(z$temperatura_bilineal_c,6))),
    rad_vecino=length(unique(round(z$radiacion_vecino_mj_m2_dia,6))),
    rad_bilineal=length(unique(round(z$radiacion_bilineal_mj_m2_dia,6))))))
resumen_unicos <- data.frame(variable=names(unicos)[3:6],
                             mediana_valores_distintos=sapply(unicos[3:6],median),
                             min=sapply(unicos[3:6],min),max=sapply(unicos[3:6],max))

# Mapas promedio y autocorrelacion descriptiva (no reemplaza el Moran residual
# que se calculara despues de definir la media en las fases 5-6).
raster_desde_espacial <- function(columna) {
  r <- mascara; values(r) <- NA_real_
  r[espacial$celda_raster] <- espacial[[columna]]; r
}
mapas <- lapply(variables,function(v) raster_desde_espacial(v)); names(mapas) <- variables
autocor <- do.call(rbind,lapply(names(mapas),function(v) {
  a <- calcular_autocorrelacion(mapas[[v]],v)
  data.frame(variable=v,n=a$n,Moran_I=a$I,Geary_C=a$C,
             Moran_esperado=-1/(a$n-1))
}))

# Controles físicos conservadores: alertan, no eliminan silenciosamente.
controles <- data.frame(
  regla=c("precipitacion < 0","temperatura fuera de [-10,45]",
          "radiacion fuera de [0,45]","altitud no finita"),
  casos=c(sum(datos$precipitacion_mm<0),
          sum(datos$temperatura_bilineal_c < -10 | datos$temperatura_bilineal_c > 45),
          sum(datos$radiacion_bilineal_mj_m2_dia < 0 |
              datos$radiacion_bilineal_mj_m2_dia > 45),
          sum(!is.finite(datos$altitud_m))))

write.csv(tabla_resumen,file.path(DIR_OUT,"estadisticas_descriptivas.csv"),row.names=FALSE)
write.csv(serie,file.path(DIR_OUT,"serie_regional_semanal.csv"),row.names=FALSE)
write.csv(climatologia,file.path(DIR_OUT,"climatologia_semanal.csv"),row.names=FALSE)
write.csv(anual,file.path(DIR_OUT,"resumen_anual.csv"),row.names=FALSE)
write.csv(espacial,file.path(DIR_OUT,"promedios_por_celda.csv"),row.names=FALSE)
write.csv(correlaciones,file.path(DIR_OUT,"correlaciones_tres_escalas.csv"),row.names=FALSE)
write.csv(comparacion_metodos,file.path(DIR_OUT,"comparacion_vecino_bilineal_eda.csv"),row.names=FALSE)
write.csv(resumen_unicos,file.path(DIR_OUT,"valores_distintos_por_metodo.csv"),row.names=FALSE)
write.csv(autocor,file.path(DIR_OUT,"autocorrelacion_espacial_descriptiva.csv"),row.names=FALSE)
write.csv(controles,file.path(DIR_OUT,"controles_fisicos.csv"),row.names=FALSE)
saveRDS(list(serie=serie,climatologia=climatologia,anual=anual,espacial=espacial,
             correlaciones=correlaciones,autocorrelacion=autocor),
        file.path(DIR_OUT,"objetos_resumen_eda.rds"))

cat("\n",strrep("=",78),"\nEDA SOBRE SOPORTE COMUN\n",sep="")
cat("Filas:",nrow(datos)," | celdas:",length(unique(datos$id_celda)),
    "| semanas:",length(unique(interaction(datos$anio,datos$semana))),"\n")
cat("\nESTADISTICAS DESCRIPTIVAS\n"); print(tabla_resumen,row.names=FALSE,digits=5)
cat("\nCONTROLES FISICOS\n"); print(controles,row.names=FALSE)
cat("\nDIFERENCIA ENTRE METODOS\n"); print(comparacion_metodos,row.names=FALSE,digits=5)
cat("\nVALORES ESPACIALES DISTINTOS POR CAPA\n"); print(resumen_unicos,row.names=FALSE)
cat("\nCORRELACION POOLED\n"); print(round(cor_pooled,3))
cat("\nCORRELACION DE PROMEDIOS ESPACIALES\n"); print(round(cor_espacial,3))
cat("\nCORRELACION DE ANOMALIAS SEMANALES\n"); print(round(cor_anomalias,3))

tema <- theme_minimal(base_size=11)
set.seed(2026); muestra <- datos[sample.int(nrow(datos),min(30000,nrow(datos))),]

dibujar_eda <- function() {
  p1 <- ggplot(datos,aes(precipitacion_mm))+geom_histogram(bins=60,fill="steelblue",color="white")+
    coord_cartesian(xlim=c(0,quantile(datos$precipitacion_mm,.99)))+
    labs(title="Precipitación semanal (hasta P99)",x="mm/semana",y="Frecuencia")+tema
  p2 <- ggplot(datos,aes(y=precipitacion_mm))+geom_boxplot(fill="steelblue",outlier.alpha=.05)+
    coord_cartesian(ylim=c(0,quantile(datos$precipitacion_mm,.99)))+
    labs(title="Precipitación: asimetría y extremos",x=NULL,y="mm/semana")+tema
  p3 <- ggplot(datos,aes(temperatura_bilineal_c))+geom_histogram(bins=50,fill="firebrick",color="white")+
    labs(title="Temperatura bilineal",x="°C",y="Frecuencia")+tema
  p4 <- ggplot(espacial,aes(altitud_m))+geom_histogram(bins=35,fill="darkgreen",color="white")+
    labs(title="Altitud de las 686 celdas",x="m",y="Frecuencia")+tema
  grid.arrange(p1,p2,p3,p4,ncol=2)

  p5 <- ggplot(serie,aes(fecha,precipitacion_mm))+geom_line(color="steelblue",linewidth=.35)+
    labs(title="Precipitación media regional 2010–2025",x=NULL,y="mm/semana")+tema
  p6 <- ggplot(serie,aes(fecha,temperatura_bilineal_c))+geom_line(color="firebrick",linewidth=.35)+
    labs(title="Temperatura media regional",x=NULL,y="°C")+tema
  p7 <- ggplot(serie,aes(fecha,radiacion_vecino_mj_m2_dia))+geom_line(color="orange3",linewidth=.35)+
    labs(title="Radiación media regional — vecino",x=NULL,y="MJ/m²/día")+tema
  p8 <- ggplot(anual,aes(anio,precipitacion_mm))+geom_line()+geom_point(color="steelblue")+
    labs(title="Media anual de precipitación semanal",x="Año",y="mm/semana")+tema
  grid.arrange(p5,p6,p7,p8,ncol=2)

  p9 <- ggplot(climatologia,aes(semana,precipitacion_mm))+geom_line(color="steelblue",linewidth=1)+
    labs(title="Ciclo climatológico de precipitación",x="Semana",y="mm/semana")+tema
  p10 <- ggplot(climatologia,aes(semana,temperatura_bilineal_c))+geom_line(color="firebrick",linewidth=1)+
    labs(title="Ciclo climatológico de temperatura",x="Semana",y="°C")+tema
  p11 <- ggplot(climatologia,aes(semana,radiacion_vecino_mj_m2_dia))+geom_line(color="orange3",linewidth=1)+
    labs(title="Ciclo climatológico de radiación",x="Semana",y="MJ/m²/día")+tema
  p12 <- ggplot(datos,aes(factor(semana),precipitacion_mm))+geom_boxplot(outlier.shape=NA,fill="lightblue")+
    coord_cartesian(ylim=c(0,quantile(datos$precipitacion_mm,.95)))+
    scale_x_discrete(breaks=as.character(seq(1,52,5)))+
    labs(title="Distribución de precipitación por semana",x="Semana",y="mm")+tema
  grid.arrange(p9,p10,p11,p12,ncol=2)

  par(mfrow=c(2,3),mar=c(2.5,2.5,3,4))
  for (v in variables) { plot(mapas[[v]],main=gsub("_"," ",v),axes=FALSE); lines(valle) }
  par(mfrow=c(1,1))

  p13 <- ggplot(muestra,aes(altitud_m,precipitacion_mm))+geom_point(alpha=.12,size=.5,color="steelblue")+
    geom_smooth(method="lm",se=FALSE,color="black")+
    coord_cartesian(ylim=c(0,quantile(muestra$precipitacion_mm,.99)))+
    labs(title="Precipitación vs altitud (muestra)",x="m",y="mm/semana")+tema
  p14 <- ggplot(muestra,aes(temperatura_bilineal_c,precipitacion_mm))+geom_point(alpha=.12,size=.5,color="firebrick")+
    geom_smooth(method="lm",se=FALSE,color="black")+
    coord_cartesian(ylim=c(0,quantile(muestra$precipitacion_mm,.99)))+
    labs(title="Precipitación vs temperatura",x="°C",y="mm/semana")+tema
  p15 <- ggplot(muestra,aes(radiacion_vecino_mj_m2_dia,precipitacion_mm))+geom_point(alpha=.12,size=.5,color="orange3")+
    geom_smooth(method="lm",se=FALSE,color="black")+
    coord_cartesian(ylim=c(0,quantile(muestra$precipitacion_mm,.99)))+
    labs(title="Precipitación vs radiación",x="MJ/m²/día",y="mm/semana")+tema
  p16 <- ggplot(muestra,aes(temperatura_vecino_c,temperatura_bilineal_c))+
    geom_point(alpha=.03,size=.25)+geom_abline(linetype=2,color="red")+
    labs(title="Temperatura: vecino vs bilineal",x="Vecino",y="Bilineal")+tema
  grid.arrange(p13,p14,p15,p16,ncol=2)

  csel <- cor_espacial
  zl <- matriz_a_largo(csel,"espacial")
  p17 <- ggplot(zl,aes(variable_1,variable_2,fill=correlacion))+
    geom_tile(color="white")+geom_text(aes(label=sprintf("%.2f",correlacion)),size=3)+
    scale_fill_gradient2(low="steelblue",mid="white",high="firebrick",midpoint=0,limits=c(-1,1))+
    labs(title="Correlación de promedios espaciales",x=NULL,y=NULL)+
    theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1))
  p18 <- ggplot(muestra,aes(radiacion_vecino_mj_m2_dia,radiacion_bilineal_mj_m2_dia))+
    geom_point(alpha=.03,size=.25)+geom_abline(linetype=2,color="red")+
    labs(title="Radiación: vecino vs bilineal",x="Vecino",y="Bilineal")+tema
  grid.arrange(p17,p18,ncol=2,widths=c(1.25,1))
}

# Visible en el historial del panel Plots.
dibujar_eda()
pdf(file.path(DIR_OUT,"EDA_completo_soporte_comun.pdf"),width=11,height=8.5,onefile=TRUE)
dibujar_eda(); dev.off()

cat("\nResultados guardados en:",DIR_OUT,"\n")
cat("El Moran mostrado es descriptivo de superficies promedio; la prueba clave\n")
cat("se hará sobre residuales después de definir el modelo de media.\n")
