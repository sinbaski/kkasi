                                        # Program to estimate in fit best model in a lognormal-Poisson regression. Both pseudo and full mle.

BIC=sigma2=rep(0,7)
sel=beta=matrix(0,7,4)
sel[1,]=c(1,2,0,0)
sel[2,]=c(1,0,3,0)
sel[3,]=c(1,0,0,4)
sel[4,]=c(1,2,3,0)
sel[5,]=c(1,0,3,4)
sel[6,]=c(1,2,0,4)
sel[7,]=c(1,2,3,4)

a=as.matrix(bulkclaim)
xx=log(a[,1])
dedx=log(a[,2])
ztx=a[,3:5]
nx=dim(ztx)[1]
ztx=cbind(rep(1,nx),ztx)

fnl=function(y){
    b=y[1:p]
    t=y[p+1]
    fnl=nx*log(t)/2+sum((xx-z%*%b)^2)/(2*t)+sum(log(pnorm((z%*%b-dedx)/sqrt(t))))
}

for(i in 1:7){
    utv=sel[i,]
    utv=utv[utv>0]
    p=length(utv)
    z=ztx[,utv]
    v=solve(t(z)%*%z)
    bet0=v%*%(t(z)%*%xx)
    sig02=sum((xx-z%*%bet0)^2)/n

    init=c(bet0,sig02)
    f=optim(init,fnl,method="BFGS",hessian=T)
    par=f$par
    beta[i,utv]=par[1:p]
    sigma2[i]=par[p+1]
    BIC[i]=2*(f$val+n*log(2*pi)/2+sum(xx))+(p+1)*log(n)
}

posx=order(BIC)[1]
utv=sel[posx,]
utvx=utv[utv>0]
bet=beta[posx,utvx]    # regression parameters for best fit
sig2=sigma2[posx]      # the sigma^2 for best fit

a=as.matrix(bulkpolicy)   # starting claim number fitting
xn=a[,6]
dedn=log(a[,2])
ztn=a[,3:5]
q0=a[,1]
nn=dim(ztn)[1]
ztn=cbind(rep(1,nn),ztn)
BIC2=rep(0,7)
gamma=matrix(0,7,4)
q=q0*pnorm((ztn[,utv]%*%bet-dedn)/sqrt(sig2))

for(i in 1:7){
    utv=sel[i,]
    utv=utv[utv>0]
    p=length(utv)
    f=glm(xn~ztn[,utv]-1,family=poisson,offset=log(q))
    gamma[i,utv]=f$coeff
    BIC2[i]=f$aic-2*length(utv)+length(utv)*log(nn)
}

posn=order(BIC2)[1]
utv=sel[posn,]
utvn=utv[utv>0]
gam=gamma[posn,utvn]    # regression parameters for best fit

zx=apply(zt[,utvx],2,mean)
zn=apply(zt[,utvn],2,mean)
d=c(0,100000,1000000)

valx=val=rep(0,3)
mux=sum(zx*bet)# Estimated risk premium
valx[1]=exp(mux+sig2/2)
valx[2]=valx[1]*pnorm((mux-log(d[2]))/sqrt(sig2)+sqrt(sig2))-d[2]*pnorm((mux-log(d[2]))/sqrt(sig2))
valx[3]=valx[1]*pnorm((mux-log(d[3]))/sqrt(sig2)+sqrt(sig2))-d[3]*pnorm((mux-log(d[3]))/sqrt(sig2))
mun=exp(sum(zn*gam))
val=round(mun*valx,0)        # Estimated risk premium

fn=function(y){                     # Starts full MLE
    b1=y[1:p1]
    t=y[p1+1]
    b2=y[(p1+2):(p1+p2+1)]
    a=ztn[,utvn]%*%b2
    fn=nx*log(t)/2+sum((xx-ztx[,utvx]%*%b1)^2)/(2*t)-sum(xn*a)+sum(q0*pnorm((ztn[,utvx]%*%b1-dedn)/sqrt(t))*exp(a))
}

p1=length(utvx)
p2=length(utvn)
init=c(bet,sig2,gam)
v=optim(init,fn,method="BFGS",hessian=T)

par=v$par
betmle=par[1:p1]         # Estimated lognormal parameters
sig2mle=par[p1+1]        # Estimated lognormal sigma^2
gammle=par[(p1+2):(p1+p2+1)]    # Estimated Poisson parameters

valx=valmle=rep(0,3)
mux=sum(zx*betmle)
valx[1]=exp(mux+sig2mle/2)
valx[2]=valx[1]*pnorm((mux-log(d[2]))/sqrt(sig2mle)+sqrt(sig2mle))-d[2]*pnorm((mux-log(d[2]))/sqrt(sig2mle))
valx[3]=valx[1]*pnorm((mux-log(d[3]))/sqrt(sig2mle)+sqrt(sig2mle))-d[3]*pnorm((mux-log(d[3]))/sqrt(sig2mle))
mun=exp(sum(zn*gammle))
valmle=round(mun*valx,0)   # Estimated risk premium



















