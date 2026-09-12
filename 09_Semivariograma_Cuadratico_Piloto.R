# =============================================================================
# 09_Semivariograma_Cuadratico_Piloto.R
#
# Ajuste de semivariograma a los residuales de la tendencia cuadrática del
# corte 2014-banda 38. Compara exponencial, esférico y gaussiano.
# =============================================================================
suppressPackageStartupMessages({library(sp); library(gstat)})
buscar_base <- function() {
  cc <- unique(normalizePath(c(getwd(),file.path(getwd(),".."),
    file.path(getwd(),"../.."),file.path(getwd(),"../../..")),mustWork=FALSE))
  ok <- vapply(cc,function(p) file.exists(file.path(p,
    "avance_team_proyecto1","reinicio_desde_cero","tendencia_flexible_resultados",
    "residuales_cuadratica_ambiental.csv")),logical(1))
  if(!any(ok)) stop("No se encontró residuales_cuadratica_ambiental.csv.")
  cc[which(ok)[1]]
}
BASE <- buscar_base()
IN <- file.path(BASE,"avance_team_proyecto1","reinicio_desde_cero",
                "tendencia_flexible_resultados","residuales_cuadratica_ambiental.csv")
OUT <- file.path(BASE,"avance_team_proyecto1","reinicio_desde_cero",
                 "semivariograma_cuadratico_resultados")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(IN)
lat0 <- mean(d$y)
d$x_km <- (d$x-mean(d$x))*111.32*cos(lat0*pi/180)
d$y_km <- (d$y-lat0)*111.32
spdf <- SpatialPointsDataFrame(coords=d[,c("x_km","y_km")],
  data=data.frame(residual=d$residual))
dist_max <- max(spDists(spdf)); cutoff <- .85*dist_max
vg <- variogram(residual~1,spdf,cutoff=cutoff,width=cutoff/15)
vg <- vg[is.finite(vg$gamma)&vg$np>0,]
write.csv(vg[,c("dist","gamma","np")],
  file.path(OUT,"semivariograma_empirico_residuales_cuadratico.csv"),row.names=FALSE)
vr <- var(d$residual); ini <- vgm(.9*vr,"Sph",cutoff/3,.1*vr)
tipos <- c(exponencial="Exp",esferico="Sph",gaussiano="Gau")
aj <- lapply(tipos,function(t) tryCatch(fit.variogram(vg,
  vgm(.9*vr,t,cutoff/3,.1*vr),fit.method=7),error=function(e) NULL))
names(aj) <- names(tipos)
tab <- do.call(rbind,lapply(names(aj),function(nm){
  a<-aj[[nm]]
  if(is.null(a)) return(data.frame(modelo=nm,nugget=NA,psill=NA,range_km=NA,SSErr=NA))
  data.frame(modelo=nm,nugget=a$psill[a$model=="Nug"][1],
    psill=a$psill[a$model!="Nug"][1],range_km=a$range[a$model!="Nug"][1],
    SSErr=attr(a,"SSErr"))
}))
tab <- tab[order(tab$SSErr),]
write.csv(tab,file.path(OUT,"comparacion_modelos_variograma_cuadratico.csv"),row.names=FALSE)
print(tab,row.names=FALSE)
pdf(file.path(OUT,"semivariograma_cuadratico_comparacion.pdf"),width=11,height=8.5)
par(mfrow=c(1,2),mar=c(4.5,4.5,3.5,1.5))
set.seed(390); pr<-combn(nrow(d),2); ii<-sample(seq_len(ncol(pr)),min(12000,ncol(pr)))
hc<-sqrt((d$x_km[pr[1,ii]]-d$x_km[pr[2,ii]])^2+(d$y_km[pr[1,ii]]-d$y_km[pr[2,ii]])^2)
gc<-(d$residual[pr[1,ii]]-d$residual[pr[2,ii]])^2/2
plot(hc,gc,pch=16,cex=.35,col=adjustcolor("grey30",.25),
  xlab="Distancia entre celdas (km)",ylab="Semivarianza",
  main="Nube: residuales cuadráticos")
points(vg$dist,vg$gamma,pch=19,col="navy")
plot(vg$dist,vg$gamma,pch=19,col="navy",xlab="Distancia (km)",
  ylab="Semivarianza residual",main="Empírico y modelos teóricos")
cc<-c(exponencial="firebrick",esferico="darkgreen",gaussiano="darkorange")
hh<-seq(0,max(vg$dist),length.out=200)
for(nm in names(aj)) if(!is.null(aj[[nm]]))
 lines(hh,variogramLine(aj[[nm]],dist_vector=hh)$gamma,col=cc[nm],lwd=2)
legend("bottomright",legend=names(cc),col=cc,lwd=2,bty="n")
dev.off()
cat("Resultados guardados en:",OUT,"\n")

