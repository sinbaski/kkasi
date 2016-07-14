objf.garch <- function(vartheta, eps,n,sig2init,petit=sqrt(.Machine$double.eps),r0=10){
  omega <- vartheta[1]
  alpha <- vartheta[2]
  beta <- vartheta[3]
  sig2<-rep(0,n)
  sig2[1]<-sig2init
  for(t in 2:n){
    sig2[t]<-omega+alpha*eps[t-1]^2+beta*sig2[t-1]
  }
  qml <- mean(eps[(r0+1):n]^2/sig2[(r0+1):n]+log(sig2[(r0+1):n]))
  qml }
#
VarAsymp<- function(omega,alpha,beta,eps,sig2init,petit,r0=10){
  n <- length(eps)
  dersigma2<-matrix(0,nrow=3,ncol=n)
  sig2<-rep(0,n)
  sig2[1]<-sig2init
  for(t in 2:n){
    vec<-c(1,eps[t-1]^2,sig2[t-1])
    sig2[t]<-omega+beta*sig2[t-1]+alpha*eps[t-1]^2
    dersigma2[1:3,t]<-vec/sig2[t]+beta*dersigma2[1:3,(t-1)]
  }
  eta <- eps[(r0+1):n]/sqrt(sig2)[(r0+1):n]
  eta <- eta/sd(eta)

  J<-dersigma2[1:3,(r0+1):n]%*%t(dersigma2[1:3,(r0+1):n])/(n-r0)
  kappa4<-mean(eta^4)

{if(kappa(J)<1/petit) inv<-solve(J) else inv<-matrix(0,nrow=3,ncol=3)}
var<-(kappa4-1)*inv
list(var=var,residus=eta)
}


estimGARCH<- function(omega,alpha,beta,eps,petit=sqrt(.Machine$double.eps),r0=10)
{
  valinit<-c(omega,alpha,beta)
  n <- length(eps)
  sig2init<-var(eps[1:min(n,5)])
  res <- nlminb(valinit,objf.garch,lower=c(petit,0,0),
                upper=c(Inf,1,1), eps=eps,n=n,sig2init=sig2init)
  omega <- res$par[1]
  alpha<- res$par[2]
  beta <- res$par[3]
  var<-VarAsymp(omega,alpha,beta,eps,sig2init,petit=sqrt(.Machine$double.eps),r0=10)
  list(coef=c(omega,alpha,beta),residus=var$residus,var=var$var)
}
omega.init<- 0
alpha.init<-0.01
beta.init<- 0
factor<-1
fitgarch<-estimGARCH(omega.init,alpha.init,beta.init,eps)
par<-fitgarch$coef
res<-fitgarch$residus
qqnorm(res)
acf(res^2)
