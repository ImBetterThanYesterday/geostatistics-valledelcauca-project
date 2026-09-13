# =============================================================================
# 16_Mapas_Top6_Modelos_Por_Corte.R
# Genera una lamina espacial por corte con las predicciones de los 6 mejores
# modelos MCO validos. Complementa la tabla de los 48 candidatos.
# No usa kriging: aqui se ve visualmente la superficie de CADA modelo de media.
# =============================================================================
buscar_raiz <- function(){
  cc <- unique(normalizePath(c(getwd(),file.path(getwd(),".."),file.path(getwd(),"../..")),mustWork=FALSE))
  ok <- vapply(cc,function(p) file.exists(file.path(p,"avance_team_proyecto1","reinicio_desde_cero","dataset_eda","dataset_final_7_variables.rds")),logical(1))
  if(!any(ok)) stop("Abra la raiz de Proyecto1 en RStudio.")
  cc[which(ok)[1]]
}
TODOS_LOS_CORTES <- TRUE
CORTE_EJEMPLO <- c(2014L,38L)
RAIZ <- buscar_raiz()
BASE <- file.path(RAIZ,"avance_team_proyecto1","reinicio_desde_cero")
OUT <- file.path(BASE,"resultados_baseR_c2","reportes_por_corte")
DATOS <- readRDS(file.path(BASE,"dataset_eda","dataset_final_7_variables.rds"))
CAND <- read.csv(file.path(BASE,"resultados_baseR_c2","modelos_candidatos.csv"),check.names=FALSE)
DEC <- read.csv(file.path(BASE,"resultados_baseR_c2","resumen_decisiones.csv"),check.names=FALSE)

inversa <- function(z,escala) if(escala=="original") z else if(escala=="raiz") pmax(z,0)^2 else pmax(expm1(z),0)
color_mapa <- function(z){ n <- 70; br <- seq(min(z,na.rm=TRUE),max(z,na.rm=TRUE),length.out=n+1); hcl.colors(n,"Blues 3",rev=TRUE)[pmax(1,pmin(n,findInterval(z,br,all.inside=TRUE)))] }

lamina_corte <- function(a,s){
  d <- DATOS[DATOS$anio==a & DATOS$semana==s,]
  co <- CAND[CAND$anio==a & CAND$semana==s & !CAND$colineal,]
  if(!nrow(d)||!nrow(co)) return(NULL)
  co <- co[order(co$RMSE),][seq_len(min(6,nrow(co))),]
  lat <- mean(d$y); d$x_km <- (d$x-mean(d$x))*111.32*cos(lat*pi/180); d$y_km <- (d$y-mean(d$y))*111.32
  carpeta <- file.path(OUT,sprintf("%d_semana_%02d",a,s)); dir.create(carpeta,recursive=TRUE,showWarnings=FALSE)
  archivo <- file.path(carpeta,sprintf("04_mapas_top6_MCO_%d_semana_%02d.png",a,s))
  png(archivo,width=2100,height=1500,res=160)
  on.exit(dev.off(),add=TRUE)
  par(mfrow=c(2,3),mar=c(3.8,3.8,3.5,.8),oma=c(0,0,2,0))
  for(i in seq_len(nrow(co))){
    d$.respuesta <- switch(co$escala[i],original=d$precipitacion_mm,raiz=sqrt(d$precipitacion_mm),log1p=log1p(d$precipitacion_mm))
    fit <- lm(as.formula(co$formula[i]),data=d)
    pred <- inversa(predict(fit,newdata=d),co$escala[i])
    plot(d$x,d$y,pch=15,cex=1.15,col=color_mapa(pred),asp=1,xlab="Longitud",ylab="Latitud",
         main=sprintf("#%d %s | %s\nRMSE = %.2f mm",i,co$modelo[i],co$escala[i],co$RMSE[i]))
  }
  mtext(sprintf("Seis mejores modelos MCO: %d, semana ISO %02d",a,s),outer=TRUE,side=3,line=.3,font=2)
  data.frame(anio=a,semana=s,archivo=archivo)
}

cortes <- unique(CAND[,c("anio","semana")])
cortes <- cortes[order(cortes$semana,cortes$anio),]
if(!TODOS_LOS_CORTES) cortes <- data.frame(anio=CORTE_EJEMPLO[1],semana=CORTE_EJEMPLO[2])
idx <- vector("list",nrow(cortes))
for(i in seq_len(nrow(cortes))){ cat(sprintf("[%d/%d] %d semana %02d\n",i,nrow(cortes),cortes$anio[i],cortes$semana[i])); idx[[i]] <- lamina_corte(cortes$anio[i],cortes$semana[i]) }
write.csv(do.call(rbind,idx),file.path(OUT,"INDICE_MAPAS_TOP6_MCO.csv"),row.names=FALSE)
cat("Listo. Las laminas estan dentro de reportes_por_corte.\n")
