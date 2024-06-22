library(WRS2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(viridis)
library(ggpubr)

win_izlase<-function(izlase,alpha){
n<-length(izlase)
r<-floor(n*alpha)+1
s<-n-floor(n*alpha)

izlase.ord<-sort(izlase)
izl_win<-c(rep(izlase.ord[r],r-1),izlase.ord[r:s],rep(izlase.ord[s],r-1))
izl_win
}
J_fun<-function(u,alpha,gamma){
w<-0
if(u<alpha){w<-0}
else if(u<=gamma){w<-(u-alpha)/(gamma-alpha)}
else if(u<=1-gamma){w<-1}
else if(u<=1-alpha){w<-(1-alpha-u)/(gamma-alpha)}
else{w<-0}
}

ST_mean<-function(izlase,alpha,gamma){
n<-length(izlase)
weights<-sapply(seq(1:n)/(n+1),J_fun,alpha,gamma)
m<-sum(weights!=0)
w.norm<-weights*m/sum(weights)
stmean<-sum(sort(izlase)*w.norm)/m
stmean
}
T_mean<-function(dati,alpha){
mean(dati,alpha)
}
trimvar<-function(izlase,alpha){
n<-length(izlase)
win.izlase<-win_izlase(izlase,alpha=alpha)
trim.sigma<-1/(n-1)/n/(1-2*alpha)^2*sum((win.izlase-mean(win.izlase))^2)
trim.sigma
}
stmeanvar<-function(dati,alpha,gamma){
dati_sort<-sort(dati)
n<-length(dati)
r<-floor(alpha*n)
m<-floor(gamma*n)
weights<-sapply(seq(1:n)/(n+1),J_fun,alpha,gamma)
const<-n/sum(weights)

if(alpha>=gamma)NA
else{

rez<-1/(m-r)*((m+m*r/n-r-m^2/n)*dati_sort[m+1]-(1+r/n)*sum(
dati_sort[(r+1):m])+
2/n*sum(dati_sort[(r+1):m]*((r+1):m))+(2-r/n)*sum(dati_sort[(n-m+1):(n-r)])+
(r*m-m^2)/n*dati_sort[n-m]-2/n*sum(dati_sort[(n-m
+1):(n-r)]*((n-m+1):(n-r))))+
1/n*(sum(dati_sort[(m+1):(n-m)])+m*dati_sort[n-m]+(m-n)*dati_sort[m+1])

infl<-c()

for(i in 1:length(dati)){
x<-dati_sort[i]
if(dati_sort[i]<=dati_sort[r]){
infl[i]<-rez
}else if((dati_sort[i]>dati_sort[r])&(dati_sort[i]<=dati_sort[m])){
infl[i]<-1/(m-r)*((i-r)*x-sum(dati_sort[(r+1):i]))-rez
}else if((dati_sort[i]>dati_sort[m])&(dati_sort[i]<=dati_sort[n-m])){
infl[i]<-1/(m-r)*((m-r)*dati_sort[m+1]-sum(dati_sort[(r+1):m]))+x-dati_sort[m+1]-rez

}else if((dati_sort[i]>dati_sort[n-m])&(dati_sort[i]<=dati_sort[n-r])){
infl[i]<-1/(m-r)*((m-r)*dati_sort[m+1]-sum(dati_sort[(r+1):m]))
+dati_sort[n-m]-dati_sort[m+1]+
1/(m-r)*((n-r)*x+(r-m)*dati_sort[n-m]-i*x+sum(dati_sort[(n-m+1):i]))-rez
}else if(dati_sort[i]>dati_sort[n-r]){
infl[i]<-1/(m-r)*((m-r)*dati_sort[m+1]-sum(dati_sort[(r+1):m]))
+dati_sort[n-m]-dati_sort[m+1]+
1/(m-r)*((r-m)*dati_sort[n-m]+sum(dati_sort[(n-m+1):(n-r)]))-
rez
}
}
sum(infl^2)*const^2/n^2
}
}


norm.appr.stmean<-function(dati,alpha,gamma,var.stmean){

stmean<-ST_mean(dati,alpha,gamma)
#var.stmean<-stmeanvar(dati,alpha,gamma)

lb<-stmean-qnorm(1-(1-0.95)/2)*sqrt(var.stmean)
ub<-stmean+qnorm(1-(1-0.95)/2)*sqrt(var.stmean)

int_length<-ub-lb
indic<-prod(c(ub,lb)-0)#atnemistovidejovertibu

c(indic,int_length)
#int_length
}
norm.appr.tmean<-function(dati,alpha,var.tmean){

tmean<-mean(dati,alpha)
#var.tmean<-trimvar(dati,alpha)

lb<-tmean-qnorm(1-(1-0.95)/2)*sqrt(var.tmean)
ub<-tmean+qnorm(1-(1-0.95)/2)*sqrt(var.tmean)

int_length<-ub-lb
indic<-prod(c(ub,lb)-0)#atnemistovidejovertibu

c(indic,int_length)
#int_length
}

emp.lik.stmean <- function(mu,dati,alpha,gamma,var.stmean){

izlase<-sort(dati)
n<-length(dati)
r<-floor(n*alpha)
weights<-sapply(seq(1:n)/(n+1),J_fun,alpha,gamma)
weights<-weights/sum(weights)
m<-n-2*r
Wni<-(izlase-mu)[(r+1):(n-r)]

lam.fun<-function(lambda)
{
sum(weights[(r+1):(n-r)]*Wni/(1+lambda*Wni))
}

lam.fun<-Vectorize(lam.fun)

gal1<-1/(-max(Wni))+0.0001
gal2<-1/(-min(Wni))-0.0001

if(gal1>=gal2){
stat<-100
}else{
sakne<-uniroot(function(lam)lam.fun(lam),c(gal1,gal2))$root
sigma1sq<-sum((weights*(izlase-ST_mean(izlase,alpha,gamma))^2))
#sigma2sq<-stmeanvar(izlase,alpha,gamma)*n
sigma2sq<-var.stmean*n
a<-sigma1sq/sigma2sq/(1-2*alpha)

stat<-2*a*sum(weights[(r+1):(n-r)]*m*log(1+sakne*Wni))}
stat
}
emp.lik.tmean<-function(mu,dati,alpha,var.tmean){

izlase<-sort(dati)
n<-length(dati)
r<-floor(n*alpha)
Wni<-(izlase-mu)[(r+1):(n-r)]
m<-n-2*r

lam.fun<-function(lambda){
1/m*sum(Wni/(1+lambda*Wni))
}

lam.fun<-Vectorize(lam.fun)

gal1<-1/(-max(Wni))+0.0001
gal2<-1/(-min(Wni))-0.0001

if(gal1>=gal2){
stat<-100
}else{

sigma1sq<-1/m*sum((izlase[(r+1):(n-r)]-mean(dati,alpha))^2)
#sigma2sq<-trimvar(izlase,alpha)*n
sigma2sq<-var.tmean*n
a<-sigma1sq/sigma2sq/(1-2*alpha)

sakne<-uniroot(function(lam)lam.fun(lam),c(gal1,gal2))$root
stat<-2*a*sum(log(1+sakne*Wni))
}
stat
}

emp.conf.intervals.stmean<-function(step=0.01,initStep=0,level=3.84,mu,dati,alpha,gamma,var.stmean){
value<-0
step1<-step
Lbeta<-mu-initStep
while(value<level){
Lbeta<-Lbeta-step1
value<-emp.lik.stmean(Lbeta,dati,alpha,gamma,var.stmean)
}
Lbeta0<-Lbeta
Lbeta1<-Lbeta+step1

tempfun<-function(beta1){
return(level-emp.lik.stmean(beta1,dati,alpha,gamma,var.stmean))
}

if(round(abs(Lbeta0-mu),2)<=0.01){
Lbeta<-mu
}else{
temp1<-uniroot(tempfun,lower=Lbeta0,upper=Lbeta1)
Lbeta<-temp1$root
}

value<-0
Ubeta<-mu+initStep
while(value<level){
Ubeta<-Ubeta+step
value<-emp.lik.stmean(Ubeta,dati,alpha,gamma,var.stmean)
}
Ubeta0<-Ubeta
Ubeta1<-Ubeta-step

if(round(abs(Ubeta0-mu),2)<=0.01){
Ubeta<-mu
}else{
temp2<-uniroot(tempfun,lower=Ubeta1,upper=Ubeta0)
Ubeta<-temp2$root
}

int_length<-Ubeta-Lbeta
indic<-prod(c(Ubeta,Lbeta)-0)
c(indic,int_length)
#indic
#int_length
}
emp.conf.intervals.tmean<-function(step=0.01,initStep=0,level=3.84,mu,dati,alpha,var.tmean){
value<-0
step1<-step
Lbeta<-mu-initStep
while(value<level){
Lbeta<-Lbeta-step1
value<-emp.lik.tmean(Lbeta,dati,alpha,var.tmean)
}
Lbeta0<-Lbeta
Lbeta1<-Lbeta+step1

tempfun<-function(beta1){
return(level-emp.lik.tmean(beta1,dati,alpha,var.tmean))
}

if(round(abs(Lbeta0-mu),2)<=0.01){
Lbeta<-mu
}else{
temp1<-uniroot(tempfun,lower=Lbeta0,upper=Lbeta1)
Lbeta<-temp1$root
}

value<-0
Ubeta<-mu+initStep
while(value<level){
Ubeta<-Ubeta+step
value<-emp.lik.tmean(Ubeta,dati,alpha,var.tmean)
}
Ubeta0<-Ubeta
Ubeta1<-Ubeta-step

if(round(abs(Ubeta0-mu),2)<=0.01){
Ubeta<-mu
}else{
temp2<-uniroot(tempfun,lower=Ubeta1,upper=Ubeta0)
Ubeta<-temp2$root
}

int_length<-Ubeta-Lbeta
indic<-prod(c(Ubeta,Lbeta)-0)
c(indic,int_length)
#indic
#int_length
}

##################Sadalijumi##############################

izlase_norm2<-function(n,prop,mu1,mu2,sd1,sd2){
y0<-rnorm(n,mu1,sd1)
y1<-rnorm(n,mu2,sd2)
flag<-rbinom(n,size=1,prob=prop)
y<-y0*flag+y1*(1-flag)
y
}
izlase_norm3<-function(n,prop1,prop2,prop3,mu1,mu2,mu3,sd1,sd2,
sd3){
mu<-c(mu1,mu2,mu3)
sd<-c(sd1,sd2,sd3)
components<-sample(1:3,prob=c(prop1,prop2,prop3),size=n,replace=
TRUE)
izlase<-rnorm(n=n,mean=mu[components],sd=sd[components])
izlase
}
izlase_norm_unif<-function(n,cpct,mu,sd,min,max){
y0<-rnorm(n,mu,sd)
y1<-runif(n,min,max)
flag<-rbinom(n,size=1,prob=cpct)
y<-y0*flag+y1*(1-flag)
y
}

##############Parametri###############################
alpha.val<-c(0.05,0.1,0.15,0.2)
gamma.val<-c(0.1,0.2,0.3,0.4)
gamma.vec<-c(gamma.val,rep(gamma.val[2:4],2),gamma.val[3:4])
alpha.vec<-c(rep(0.05,4),rep(0.1,3),rep(0.15,3),rep(0.2,2))
lines<-c(rep(1,4),2)
gamma.indic<-c(0,3,6,8)
cov.acc.tab<-data.frame()
int.len.tab<-data.frame()
#var.tab<-data.frame()
########################################################
N<-10000
# n <- 50
for(n in c(50,100,200,500)){

TM.norm.indic<-matrix(data=NA,N,4)
TM.norm.len<-matrix(data=NA,N,4)
TM.emp.indic<-matrix(data=NA,N,4)
TM.emp.len<-matrix(data=NA,N,4)

STM.norm.indic<-matrix(data=NA,N,12)
STM.norm.len<-matrix(data=NA,N,12)
STM.emp.indic<-matrix(data=NA,N,12)
STM.emp.len<-matrix(data=NA,N,12)

#TM.var<-matrix(data=NA,N,4)
#STM.var<-matrix(data=NA,N,12)
# dati <- rnorm(100)
  for(i in 1:N){
  dati<-izlase_norm2(n,0.8,0,0,1,25)

  for(a in alpha.val){
  Tmu<-mean(dati,trim=a)
  var.tmean<-trimvar(dati,alpha=a)
  norm.tmean<-norm.appr.tmean(dati,alpha=a,var.tmean=var.tmean)
  emp.tmean<-emp.conf.intervals.tmean(mu=Tmu,dati=dati,alpha=a,var.tmean=var.tmean)
  j<-match(a,alpha.val)

  TM.norm.indic[i,j]<-norm.tmean[1]
  TM.norm.len[i,j]<-norm.tmean[2]
  TM.emp.indic[i,j]<-emp.tmean[1]
  TM.emp.len[i,j]<-emp.tmean[2]
#TM.var[i,j]<-var.tmean

    for(g in gamma.val){
      if(a>=g|abs(g-a)<0.01)
        {
        next
        }
        STmu<-ST_mean(dati,a,g)
        var.stmean<-stmeanvar(dati,a,g)

        norm.stmean<-norm.appr.stmean(dati,alpha=a,gamma=g,var.stmean=var.stmean)
        emp.stmean<-emp.conf.intervals.stmean(mu=STmu,dati=dati,alpha=a,gamma=g,var.stmean=var.stmean)

        k<-match(g,gamma.val)
        l<-gamma.indic[j]

        STM.norm.indic[i,k+l]<-norm.stmean[1]
        STM.norm.len[i,k+l]<-norm.stmean[2]
        STM.emp.indic[i,k+l]<-emp.stmean[1]
        STM.emp.len[i,k+l]<-emp.stmean[2]

      }
    }
  }
  print(n)

##############Parklajumaprecizitate#######################

#TMnormalaaproksimacija,parklajumaprecizitate
  colnames(TM.norm.indic)<-paste(alpha.val)
  TM.norm.indic<-as.data.frame(TM.norm.indic)

  TM.tab.norm<-TM.norm.indic %>%
    gather(key=alpha,value=value)%>%
    group_by(alpha)%>%
    summarise(cov.acc=length(value[value<0])/N, .groups='drop')%>%
    mutate(method.gr1='tmean',method.gr2='norm')

#TMempiriskaticamiba,parklajumaprecizitate
  colnames(TM.emp.indic)<-paste(alpha.val)
  TM.emp.indic<-as.data.frame(TM.emp.indic)

  TM.tab.emp<-TM.emp.indic%>%gather(key=alpha,value=value)%>%
    group_by(alpha)%>%
    summarise(cov.acc=length(value[value<0])/N,.groups='drop')%>%
    mutate(method.gr1='tmean',
    method.gr2='emp')

#STMnormalaaproksimacija,parklajumaprecizitate
  colnames(STM.norm.indic)<-paste(alpha.vec,gamma.vec,sep="_")
  STM.norm.indic<-as.data.frame(STM.norm.indic)

  STM.tab.norm<-STM.norm.indic%>%gather(key=param,value=value)%>%
    group_by(param)%>%
    summarise(cov.acc=length(value[value<0])/N,.groups='drop')%>%
    separate(param,c("alpha","gamma"),"_")%>%
    mutate(method.gr1=paste("gamma",gamma,sep="_"),
    method.gr2='norm')

#STMempiriskaticamiba,parklajumaprecizitate
  colnames(STM.emp.indic)<-paste(alpha.vec,gamma.vec,sep="_")
  STM.emp.indic<-as.data.frame(STM.emp.indic)

  STM.tab.emp<-STM.emp.indic%>%gather(key=param,value=value)%>%
    group_by(param)%>%
    summarise(cov.acc=length(value[value<0])/N,.groups='drop')%>%
    separate(param,c("alpha","gamma"),"_")%>%
    mutate(method.gr1=paste("gamma",gamma,sep="_"),
    method.gr2='emp')

##############Ticamibasintervalugarums#######################

#TMnormalaaproksimacija
  colnames(TM.norm.len)<-paste(alpha.val)
  TM.norm.len<-as.data.frame(TM.norm.len)

  TM.tab.norm.len<-TM.norm.len%>%gather(key=alpha,value=value)%>%
    group_by(alpha)%>%
    summarise(int.len=mean(value),.groups='drop')%>%
    mutate(method.gr1='tmean',method.gr2='norm')

#TMempiriskaticamiba
  colnames(TM.emp.len)<-paste(alpha.val)
  TM.emp.len<-as.data.frame(TM.emp.len)

  TM.tab.emp.len<-TM.emp.len%>%gather(key=alpha,value=value)%>%
    group_by(alpha)%>%
    summarise(int.len=mean(value),.groups='drop')%>%
    mutate(method.gr1='tmean',
          method.gr2='emp')

#STMnormalaaproksimacija
  colnames(STM.norm.len)<-paste(alpha.vec,gamma.vec,sep="_")
  STM.norm.len<-as.data.frame(STM.norm.len)

  STM.tab.norm.len<-STM.norm.len%>%gather(key=param,value=value)%>%
    group_by(param)%>%
    summarise(int.len=mean(value),.groups='drop')%>%
    separate(param,c("alpha","gamma"),"_")%>%
    mutate(method.gr1=paste("gamma",gamma,sep="_"),
           method.gr2='norm')

#STMempiriskaticamiba
  colnames(STM.emp.len)<-paste(alpha.vec,gamma.vec,sep="_")
  STM.emp.len<-as.data.frame(STM.emp.len)

  STM.tab.emp.len<-STM.emp.len%>%gather(key=param,value=value)%>%
  group_by(param)%>%
  summarise(int.len=mean(value),.groups='drop')%>%
  separate(param,c("alpha","gamma"),"_")%>%
  mutate(method.gr1=paste("gamma",gamma,sep="_"),
  method.gr2='emp')
  
  cov.acc.tab.n<-STM.tab.norm%>%
  bind_rows(STM.tab.emp,TM.tab.norm,TM.tab.emp)%>%
  mutate(N=n)
  cov.acc.tab<-cov.acc.tab%>%bind_rows(cov.acc.tab.n)


  int.len.tab.n<-STM.tab.norm.len%>%
  bind_rows(STM.tab.emp.len,TM.tab.norm.len,TM.tab.emp.len)%>%
  mutate(N=n)
  int.len.tab<-int.len.tab%>%bind_rows(int.len.tab.n)
  plot_list_cov_acc=list()
  plot_list_int_len=list()


  for(i in 1:4){

    j<-c(50,100,200,500)

    cov.acc<-cov.acc.tab%>%
      filter(N==j[i])%>%
      mutate(alpha=as.numeric(alpha))%>%
      ggplot(aes(alpha,cov.acc))+
      geom_point(size=2,aes(linetype=method.gr2,color=method.gr1,shape=method.gr2))+
      geom_line(size=1,aes(linetype=method.gr2,color=method.gr1))+
      scale_color_manual(values=c(viridis_pal()(4),"purple"),labels=c(
      lapply(sprintf('$STM\\gamma=%#.1f$',gamma.val),TeX),"TM"))+
      geom_hline(yintercept=0.95,linetype="dashed",color="black")+
      ggtitle(paste("N=",j[i]))+
      xlab(TeX('$\\alpha$'))+
      ylab('āāprkljumaāprecizitte')+
      theme_minimal()+
      theme(plot.title=element_text(hjust=0.5))+
      labs(shape="īTicambasāintervlaveids",linetype="īTicambasāintervlaveids",color="ēēāNovrttjs")

    int.len<-int.len.tab%>%
      filter(N==j[i])%>%
      mutate(alpha=as.numeric(alpha))%>%
      ggplot(aes(alpha,int.len))+
      geom_point(size=2,aes(linetype=method.gr2,color=method.gr1,shape=method.gr2))+
      geom_line(size=1,aes(linetype=method.gr2,color=method.gr1))+
      scale_color_manual(values=c(viridis_pal()(4),"purple"),labels=c(
      lapply(sprintf('$STM\\gamma=%#.1f$',gamma.val),TeX),"TM"))+
      ggtitle(paste("N=",j[i]))+
      xlab(TeX('$\\alpha$'))+
      ylab('ēvidjaisāintervlagarums')+
      theme_minimal()+
      theme(plot.title=element_text(hjust=0.5))+
      labs(shape="īTicambasāintervlaveids",linetype="īTicambasāintervlaveids",color="ēēāNovrttjs")

    plot_list_cov_acc[[i]]=cov.acc
    plot_list_int_len[[i]]=int.len
}

ggarrange(
plot_list_cov_acc[[1]],plot_list_int_len[[1]],
plot_list_cov_acc[[2]],plot_list_int_len[[2]],
plot_list_cov_acc[[3]],plot_list_int_len[[3]],
plot_list_cov_acc[[4]],plot_list_int_len[[4]],
common.legend=TRUE,legend="bottom",ncol=2,nrow=4)
}
