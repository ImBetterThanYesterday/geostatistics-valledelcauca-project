# Implementación base R: MCO -> variograma -> Sigma -> MCG -> kriging.
# C0 son tres cortes que pertenecen al total de 832.
MODO <- "C2" # C0, C1 (2010), C2 (2011-2025)
B_MORAN <- 199L; N_LAGS <- 12L; VIF_LIM <- 10; MEJORA_MIN <- .01

raiz <- function() {
  p <- unique(normalizePath(c(getwd(),"..","../..","../../.."), mustWork=FALSE))
  ok <- vapply(p, function(z) file.exists(file.path(z,"avance_team_proyecto1","reinicio_desde_cero","dataset_eda","dataset_final_7_variables.rds")), logical(1))
  if (!any(ok)) stop("No se encontró el dataset final")
  p[which(ok)[1]]
}
BASE <- raiz()
DATOS <- readRDS(file.path(BASE,"avance_team_proyecto1","reinicio_desde_cero","dataset_eda","dataset_final_7_variables.rds"))
OUT <- file.path(BASE,"avance_team_proyecto1","reinicio_desde_cero",paste0("resultados_baseR_",tolower(MODO)))
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

RHS <- list(
 M0="1", Mxy="x_km+y_km", Malt="altitud_m", MT="temperatura_bilineal_c", MR="radiacion_bilineal_mj_m2_dia",
 MTR="temperatura_bilineal_c+radiacion_bilineal_mj_m2_dia", MTA="temperatura_bilineal_c+altitud_m", MRA="radiacion_bilineal_mj_m2_dia+altitud_m",
 MTRA="temperatura_bilineal_c+radiacion_bilineal_mj_m2_dia+altitud_m", Mxy_alt="x_km+y_km+altitud_m",
 Mxy_T="x_km+y_km+temperatura_bilineal_c", Mxy_R="x_km+y_km+radiacion_bilineal_mj_m2_dia",
 Mxy_TR="x_km+y_km+temperatura_bilineal_c+radiacion_bilineal_mj_m2_dia",
 Mxy_TA="x_km+y_km+temperatura_bilineal_c+altitud_m", Mxy_RA="x_km+y_km+radiacion_bilineal_mj_m2_dia+altitud_m",
 Mxy_TRA="x_km+y_km+temperatura_bilineal_c+radiacion_bilineal_mj_m2_dia+altitud_m")

dist_xy <- function(x,y) sqrt(outer(x,x,"-")^2+outer(y,y,"-")^2)
inv_esc <- function(z,e) if(e=="original") z else if(e=="raiz") pmax(z,0)^2 else pmax(expm1(z),0)
metricas <- function(y,p) { ok<-is.finite(y)&is.finite(p); er<-y[ok]-p[ok]; c(n=sum(ok),MAE=mean(abs(er)),RMSE=sqrt(mean(er^2)),sesgo=mean(p[ok]-y[ok]),R2=1-sum(er^2)/sum((y[ok]-mean(y[ok]))^2)) }
vif <- function(m) { X<-model.matrix(m)[,-1,drop=FALSE]; if(ncol(X)<2) return(if(ncol(X))1 else NA_real_); max(sapply(seq_len(ncol(X)),function(j){r<-summary(lm(X[,j]~X[,-j,drop=FALSE]))$r.squared;1/max(1-r,.Machine$double.eps)})) }
loo_mco <- function(m,e) { y<-model.response(model.frame(m)); inv_esc(y-residuals(m)/pmax(1-hatvalues(m),.Machine$double.eps),e) }

moran <- function(z,D) {
 n<-length(z); W<-matrix(0,n,n); for(i in seq_len(n)) W[i,order(D[i,])[2:5]]<-1/4
 q<-z-mean(z); I<-n/sum(W)*sum(W*outer(q,q))/sum(q^2); set.seed(1000+n)
 pp<-replicate(B_MORAN,{u<-sample(q);n/sum(W)*sum(W*outer(u,u))/sum(u^2)})
 c(I=I,p=(1+sum(pp>=I))/(B_MORAN+1))
}

emp_variograma <- function(r,D) {
 ij<-which(upper.tri(D),arr.ind=TRUE); h<-D[upper.tri(D)]; g<-(r[ij[,1]]-r[ij[,2]])^2/2; co<-quantile(h,.85)
 gr<-cut(h[h<=co],breaks=seq(0,co,length.out=N_LAGS+1),include.lowest=TRUE)
 data.frame(h=as.numeric(tapply(h[h<=co],gr,mean)),gamma=as.numeric(tapply(g[h<=co],gr,mean)),np=as.integer(tapply(g[h<=co],gr,length)))
}
gamma_modelo <- function(h,t,f) { u<-h/t[3]; z<-switch(f,exponencial=1-exp(-u),esferico=ifelse(u<1,1.5*u-.5*u^3,1),gaussiano=1-exp(-u^2));t[1]+t[2]*z }
ajustar_vg <- function(vg,r) {
 vg<-vg[complete.cases(vg),]; vr<-var(r); ini<-unique(quantile(vg$h,c(.2,.5,.8))); ans<-list();k<-1
 for(f in c("exponencial","esferico","gaussiano")) for(met in c("MCO","MCP")) {
  best<-NULL; val<-Inf
  for(a in ini) for(nu in c(0,.1*vr)) for(ps in c(.5*vr,vr)) {
   ob<-function(t){e<-vg$gamma-gamma_modelo(vg$h,t,f);if(met=="MCP")sum(vg$np*e^2)else sum(e^2)}
   o<-tryCatch(optim(c(nu,ps,a),ob,method="L-BFGS-B",lower=c(0,0,min(vg$h)/10),upper=c(2*vr,3*vr,2*max(vg$h))),error=function(e)NULL)
   if(!is.null(o)&&o$value<val){best<-o;val<-o$value}
  }
  if(!is.null(best)){ans[[k]]<-data.frame(familia=f,metodo=met,nugget=best$par[1],psill=best$par[2],range_km=best$par[3],SSE=best$value);k<-k+1}
 }
 tab<-do.call(rbind,ans); sel<-tab[tab$metodo=="MCO",]; list(tabla=tab,mejor=sel[which.min(sel$SSE),])
}
sigma_vg <- function(D,v) { u<-D/v$range_km; z<-switch(as.character(v$familia),exponencial=1-exp(-u),esferico=ifelse(u<1,1.5*u-.5*u^3,1),gaussiano=1-exp(-u^2)); S<-v$psill*(1-z);diag(S)<-v$nugget+v$psill;(S+t(S))/2 }
gls <- function(X,y,S) { Q<-solve(S); V<-solve(t(X)%*%Q%*%X); b<-V%*%t(X)%*%Q%*%y;list(b=drop(b),se=sqrt(diag(V)),fit=drop(X%*%b),Q=Q) }

procesar <- function(d,a,s) {
 out<-data.frame(anio=a,semana=s,n_celdas=nrow(d),estado="INICIADO")
 req<-c("id_celda","x","y","altitud_m","precipitacion_mm","temperatura_bilineal_c","radiacion_bilineal_mj_m2_dia")
 if(nrow(d)!=686||length(unique(d$id_celda))!=686||anyNA(d[,req])||any(d$precipitacion_mm<0)){out$estado<-"FALLIDO_CONTROL";return(list(resumen=out))}
 la<-mean(d$y);d$x_km<-(d$x-mean(d$x))*111.32*cos(la*pi/180);d$y_km<-(d$y-mean(d$y))*111.32
 esc<-list(original=d$precipitacion_mm,raiz=sqrt(d$precipitacion_mm),log1p=log1p(d$precipitacion_mm)); cc<-list();k<-1
 for(e in names(esc)){d$.respuesta<-esc[[e]];for(nm in names(RHS)){m<-lm(as.formula(paste(".respuesta~",RHS[[nm]])),d);v<-vif(m);mt<-metricas(d$precipitacion_mm,loo_mco(m,e));cc[[k]]<-data.frame(anio=a,semana=s,escala=e,modelo=nm,formula=paste(deparse(formula(m)),collapse=" "),vif_max=v,colineal=is.finite(v)&&v>VIF_LIM,t(mt));k<-k+1}}
 cc<-do.call(rbind,cc); va<-cc[!cc$colineal,];if(!nrow(va)){out$estado<-"FALLIDO_COLINEALIDAD";return(list(resumen=out,candidatos=cc))}
 g<-va[which.min(va$RMSE),];d$.respuesta<-esc[[g$escala]];m<-lm(as.formula(g$formula),d);r<-residuals(m);D<-dist_xy(d$x_km,d$y_km);mo<-moran(r,D)
 out<-cbind(out,g[,c("escala","modelo","formula","vif_max","MAE","RMSE","sesgo","R2")]);out$moran_I<-mo["I"];out$moran_p<-mo["p"]
 if(mo["p"]>=.05){out$estado<-"FINAL_MCO";return(list(resumen=out,candidatos=cc))}
 av<-ajustar_vg(emp_variograma(r,D),r);fit<-av$mejor;out<-cbind(out,fit[,c("familia","nugget","psill","range_km","SSE")])
 if(fit$psill/(fit$nugget+fit$psill)<.05){out$estado<-"FINAL_MCO_NUGGET_PURO";return(list(resumen=out,candidatos=cc,variogramas=av$tabla))}
 S<-sigma_vg(D,fit);if(inherits(try(chol(S),silent=TRUE),"try-error")){out$estado<-"REVISAR_COVARIANZA";return(list(resumen=out,candidatos=cc,variogramas=av$tabla))}
 X<-model.matrix(m);z<-gls(X,d$.respuesta,S);rr<-d$.respuesta-z$fit;pr_esc<-z$fit+rr-drop(z$Q%*%rr)/diag(z$Q);pr<-inv_esc(pr_esc,g$escala);mk<-metricas(d$precipitacion_mm,pr)
 out$tipo_kriging<-if(g$modelo=="M0")"ordinario"else"universal";out$RMSE_mcg_kriging<-mk["RMSE"];out$MAE_mcg_kriging<-mk["MAE"];out$sesgo_mcg_kriging<-mk["sesgo"];out$R2_predictivo_mcg_kriging<-mk["R2"];out$n_validacion_kriging<-mk["n"];out$cobertura95_escala<-mean(abs(d$.respuesta-pr_esc)<=1.96*sqrt(pmax(1/diag(z$Q),0)));out$validacion<-"LOO_condicional_covarianza_fija";out$estado<-if(mk["RMSE"]<g$RMSE*(1-MEJORA_MIN))"FINAL_KRIGING"else"FINAL_MCG"
 co<-data.frame(anio=a,semana=s,termino=colnames(X),beta_mco=coef(m),se_mco=summary(m)$coef[,2],beta_mcg=z$b,se_mcg=z$se)
 pp<-data.frame(anio=a,semana=s,id_celda=d$id_celda,x=d$x,y=d$y,observado_mm=d$precipitacion_mm,prediccion_mm=pr,varianza_condicional_escala=1/diag(z$Q))
 list(resumen=out,candidatos=cc,variogramas=av$tabla,coeficientes=co,predicciones=pp)
}

 cortes<-switch(MODO,C0=data.frame(anio=c(2020L,2014L,2019L),semana=c(29L,38L,46L)),C1=expand.grid(anio=2010L,semana=1:52),C2=expand.grid(anio=2011:2025,semana=1:52))
res<-vector("list",nrow(cortes))
juntar<-function(n){z<-lapply(res,`[[`,n);z<-z[!vapply(z,is.null,logical(1))];if(!length(z))return(NULL);nm<-unique(unlist(lapply(z,names)));z<-lapply(z,function(a){a[setdiff(nm,names(a))]<-NA;a[,nm,drop=FALSE]});do.call(rbind,z)}
guardar_parcial<-function(){write.csv(juntar("resumen"),file.path(OUT,"checkpoint_parcial.csv"),row.names=FALSE)}
for(i in seq_len(nrow(cortes))){a<-cortes$anio[i];s<-cortes$semana[i];cat(sprintf("[%d/%d] %d semana %02d\n",i,nrow(cortes),a,s));res[[i]]<-procesar(DATOS[DATOS$anio==a&DATOS$semana==s,],a,s);if(i%%10==0)guardar_parcial()}
write.csv(juntar("resumen"),file.path(OUT,"resumen_decisiones.csv"),row.names=FALSE)
write.csv(juntar("candidatos"),file.path(OUT,"modelos_candidatos.csv"),row.names=FALSE)
if(!is.null(juntar("variogramas")))write.csv(juntar("variogramas"),file.path(OUT,"semivariogramas.csv"),row.names=FALSE)
if(!is.null(juntar("coeficientes")))write.csv(juntar("coeficientes"),file.path(OUT,"mco_mcg.csv"),row.names=FALSE)
saveRDS(juntar("predicciones"),file.path(OUT,"predicciones.rds"))
cat("Finalizado:",OUT,"\n")
