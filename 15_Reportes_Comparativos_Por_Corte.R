# =============================================================================
# Reportes comparativos por corte. Usa resultados existentes; no reentrena.
# =============================================================================
buscar_raiz <- function() {
  cc <- unique(normalizePath(c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../..")), mustWork=FALSE))
  ok <- vapply(cc, function(p) file.exists(file.path(p, "avance_team_proyecto1", "reinicio_desde_cero", "resultados_baseR_c2", "modelos_candidatos.csv")), logical(1))
  if (!any(ok)) stop("Abra la raiz de Proyecto1 en RStudio.")
  cc[which(ok)[1]]
}

# TRUE genera una carpeta por cada uno de los 780 cortes.
GENERAR_TODOS <- TRUE
CORTE_EJEMPLO <- c(2014L, 38L)
RAIZ <- buscar_raiz()
RUTA <- file.path(RAIZ, "avance_team_proyecto1", "reinicio_desde_cero", "resultados_baseR_c2")
SALIDA <- file.path(RUTA, "reportes_por_corte")
dir.create(SALIDA, recursive=TRUE, showWarnings=FALSE)
resumen <- read.csv(file.path(RUTA, "resumen_decisiones.csv"), check.names=FALSE)
candidatos <- read.csv(file.path(RUTA, "modelos_candidatos.csv"), check.names=FALSE)
predicciones <- readRDS(file.path(RUTA, "predicciones.rds"))

guardar_png <- function(archivo, codigo, ancho=1500, alto=950) {
  png(archivo, width=ancho, height=alto, res=150); on.exit(dev.off(), add=TRUE)
  eval.parent(substitute(codigo))
}
colores_mapa <- function(z, residual=FALSE) {
  n <- 80
  br <- if(residual) { q <- max(abs(z),na.rm=TRUE); seq(-q,q,length.out=n+1) } else seq(min(z,na.rm=TRUE),max(z,na.rm=TRUE),length.out=n+1)
  pal <- if(residual) hcl.colors(n,"Blue-Red 3") else hcl.colors(n,"Blues 3",rev=TRUE)
  pal[pmax(1,pmin(n,findInterval(z,br,all.inside=TRUE)))]
}

# Tabla global: los 37.440 candidatos, incluyendo los excluidos por VIF.
global <- candidatos
global$valido <- !global$colineal
global$rank_RMSE_en_corte <- ave(global$RMSE, interaction(global$anio,global$semana), FUN=function(z) rank(z,ties.method="first"))
gana <- resumen[,c("anio","semana","modelo","escala","estado")]
names(gana)[3:4] <- c("modelo_ganador","escala_ganadora")
global <- merge(global,gana,by=c("anio","semana"),all.x=TRUE)
global$ganador_modelo_medio <- global$valido & global$modelo==global$modelo_ganador & global$escala==global$escala_ganadora
write.csv(global,file.path(SALIDA,"TABLA_GLOBAL_37440_MODELOS.csv"),row.names=FALSE)

# Una grafica global: distribucion de RMSE por modelo + transformacion.
vv <- global[global$valido,]; vv$etiqueta <- paste(vv$modelo,vv$escala,sep=" | ")
orden <- names(sort(tapply(vv$RMSE,vv$etiqueta,median,na.rm=TRUE)))
vv$etiqueta <- factor(vv$etiqueta,levels=orden)
guardar_png(file.path(SALIDA,"00_comparacion_global_48_modelos.png"), {
  par(mar=c(12,5,4,2)); boxplot(RMSE~etiqueta,data=vv,las=2,cex.axis=.65,col="#74ADD1",outline=FALSE,
    main="Comparacion global de los candidatos validos",ylab="RMSE LOO (mm)")
  mtext("Cada caja resume un modelo + transformacion en todos los cortes. Menor es mejor.",side=3,line=.3,cex=.8)
},2400,1350)

reporte_corte <- function(a,s) {
  etiqueta <- sprintf("%d_semana_%02d",a,s); carpeta <- file.path(SALIDA,etiqueta)
  dir.create(carpeta,recursive=TRUE,showWarnings=FALSE)
  comp <- candidatos[candidatos$anio==a & candidatos$semana==s,]
  dec <- resumen[resumen$anio==a & resumen$semana==s,]
  if(!nrow(comp)||!nrow(dec)) return(NULL)
  comp$valido <- !comp$colineal
  comp$rank_RMSE <- NA_integer_; ii <- which(comp$valido); comp$rank_RMSE[ii] <- rank(comp$RMSE[ii],ties.method="first")
  comp$ganador <- comp$valido & comp$modelo==dec$modelo[1] & comp$escala==dec$escala[1]
  comp$estado_candidato <- ifelse(!comp$valido,"EXCLUIDO: VIF > 10",ifelse(comp$ganador,"GANADOR MCO","VALIDO"))
  comp <- comp[order(!comp$valido,comp$RMSE),]
  write.csv(comp,file.path(carpeta,paste0("tabla_modelos_",etiqueta,".csv")),row.names=FALSE)

  # Grafica de TODOS los candidatos: azul valido, rojo ganador, gris VIF >10.
  guardar_png(file.path(carpeta,paste0("01_comparacion_modelos_",etiqueta,".png")), {
    par(mar=c(12,5,5,2)); labs <- paste(comp$modelo,comp$escala,sep="\n")
    cols <- ifelse(!comp$valido,"grey80",ifelse(comp$ganador,"#D73027","#4575B4"))
    barplot(comp$RMSE,names.arg=labs,las=2,cex.names=.62,col=cols,border=NA,ylab="RMSE LOO (mm)",
      main=sprintf("Modelos probados: %d, semana ISO %02d",a,s))
    abline(h=dec$RMSE[1],lty=2,lwd=2,col="#D73027")
    legend("topright",legend=c("Valido","Ganador","Excluido: VIF > 10"),fill=c("#4575B4","#D73027","grey80"),bty="n")
  })
  top <- comp[comp$valido,][seq_len(min(6,sum(comp$valido))),]
  guardar_png(file.path(carpeta,paste0("02_top6_metricas_",etiqueta,".png")), {
    par(mfrow=c(1,3),mar=c(10,4,4,1)); labs <- paste(top$modelo,top$escala,sep="\n")
    barplot(top$RMSE,names.arg=labs,las=2,col="#4575B4",main="RMSE LOO",ylab="mm")
    barplot(top$MAE,names.arg=labs,las=2,col="#74ADD1",main="MAE LOO",ylab="mm")
    barplot(top$R2,names.arg=labs,las=2,col="#ABD9E9",main="R2 predictivo MCO",ylab="R2")
  },1800,950)
  sub <- predicciones[predicciones$anio==a & predicciones$semana==s,]
  if(nrow(sub)) {
    sub$residual <- sub$observado_mm-sub$prediccion_mm
    guardar_png(file.path(carpeta,paste0("03_resultado_kriging_",etiqueta,".png")), {
      par(mfrow=c(2,2),mar=c(4,4,3,1),oma=c(0,0,2,0))
      plot(sub$x,sub$y,pch=15,cex=1.25,col=colores_mapa(sub$observado_mm),asp=1,xlab="Longitud",ylab="Latitud",main="Precipitacion observada")
      plot(sub$x,sub$y,pch=15,cex=1.25,col=colores_mapa(sub$prediccion_mm),asp=1,xlab="Longitud",ylab="Latitud",main="Prediccion por kriging")
      plot(sub$x,sub$y,pch=15,cex=1.25,col=colores_mapa(sub$residual,TRUE),asp=1,xlab="Longitud",ylab="Latitud",main="Residual: observado - predicho")
      plot(sub$observado_mm,sub$prediccion_mm,pch=16,col=rgb(.1,.3,.7,.5),xlab="Precipitacion observada (mm)",ylab="Prediccion LOO (mm)",main="Validacion punto a punto")
      abline(0,1,lty=2,lwd=2,col="red3"); mtext(sprintf("Resultado final: %d, semana ISO %02d",a,s),outer=TRUE,side=3,line=.3,font=2)
    },1800,1500)
  }
  data.frame(anio=a,semana=s,estado=dec$estado[1],modelo_ganador=dec$modelo[1],escala=dec$escala[1],RMSE_MCO=dec$RMSE[1],RMSE_kriging=dec$RMSE_mcg_kriging[1])
}

cortes <- resumen[,c("anio","semana")]
if(!GENERAR_TODOS) cortes <- data.frame(anio=CORTE_EJEMPLO[1],semana=CORTE_EJEMPLO[2])
indice <- vector("list",nrow(cortes))
for(i in seq_len(nrow(cortes))) { cat(sprintf("[%d/%d] %d semana %02d\n",i,nrow(cortes),cortes$anio[i],cortes$semana[i])); indice[[i]] <- reporte_corte(cortes$anio[i],cortes$semana[i]) }
write.csv(do.call(rbind,indice),file.path(SALIDA,"INDICE_REPORTES_POR_CORTE.csv"),row.names=FALSE)
cat("Listo. Revise:",SALIDA,"\n")
