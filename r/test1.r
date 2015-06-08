library(RMySQL);
library(alabama);
library(forecast);
rm(list=ls());

source("libxxie.r");

my.fun <- function(coef, R) {
    auto <- acf(R %*% coef, plot = FALSE)$acf[-1];
    p <- length(auto);
    auto <- auto * exp(-(0:(p-1)/5));
    auto <- auto[which(abs(auto) > 2/sqrt(dim(R)[1]))];
    return(-sum(abs(auto)));
}

heq <- function(coef, R) {
    return(sum(abs(coef)) - 1);
}

day2 = '2015-05-28';
day1 = '2011-01-01';

database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host=Sys.getenv("PB"));

assetSet <- "indices";
results = dbSendQuery(database, sprintf("select tblname from %s;", assetSet));
tables <- fetch(results, n=-1)[[1]];
p = length(tables);
dbClearResult(results);
dbDisconnect(database);
data <- getAssetReturns(day1, day2, tables);
R = matrix(unlist(data[, -1]), nrow=dim(data)[1], byrow=FALSE);
T = dim(R)[1];

A <- rep(1/p, p);
A[1] <- 1;
result <- auglag(par=A, fn=my.fun, R=R, heq=heq);
ret <- R %*% result$par;

library(rugarch);
## For the generalized hyperbolic distribution
## which is a normal mean-variance mixture, there
## exists location and scale invariant parametrization,
## namely ((lambda,zeta), rho), where (lambda, zeta)
## constitute the shape parameters and rho is the skewness
## parameter.
##
## When lambda is fixed
## zeta = delta * (alpha + beta)^{1/2}
## rho = beta/alpha
##
## For the generalized hyperbolic skewed t distribution
## alpha 
spec <- ugarchspec(mean.model=list(armaOrder=c(2, 1)),
                   distribution.model="ghst",
                   variance.model=list(
                       model="gjrGARCH",
                       garchOrder=c(1, 1)
                   )
                   );
model <- ugarchfit(spec=spec, data=ret);
library(parallel);
cluster <- makePSOCKcluster(4);
roll <- ugarchroll(spec, ret, n.start=floor(T*0.9), refit.every=60,
                   refit.window="moving", window.size=200,
                   solver="hybrid", calculate.VaR=TRUE,
                   VaR.alpha=0.05, cluster=cluster, keep.coef=TRUE);
stopCluster(cluster);

## acov <- acf(ret, type="covariance", plot=FALSE);
## inno <- inferInnovations(ret);


## Akaike <- matrix(nrow=5,ncol=5);
## Beysian <- matrix(nrow=5,ncol=5);
## for (p in 1:5) {
##     for (q in 1:5) {
##         spec <- ugarchspec(mean.model=list(armaOrder=c(p, q)),
##                            distribution.model="ghst",
##                            variance.model=list(
##                                model="sGARCH",
##                                garchOrder=c(1, 1)
##                            )
##                            );
##         model <- ugarchfit(spec=spec, data=ret);
##         if (0 == convergence(model)) {
##             Akaike[p,q] <- infocriteria(model)[1];
##             Beysian[p,q] <- infocriteria(model)[2];
##         }
##     }
## }


