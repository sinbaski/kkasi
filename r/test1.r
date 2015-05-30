library(RMySQL);
library(alabama);
library(forecast);
rm(list=ls());

source("libxxie.r");

my.fun <- function(coef, R) {
    auto <- acf(R %*% coef, plot = FALSE)$acf[-1];
    p <- length(auto);
    au <- auto * exp(-(0:p)/5);
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
ret <- getAssetReturns(day1, day2, tables);
R = matrix(unlist(ret[, -1]), nrow=dim(ret)[1], byrow=FALSE);
T = dim(R)[1];

A <- rep(1/p, p);
A[1] <- 1;
result <- auglag(par=A, fn=my.fun, R=R, heq=heq);
ret <- R %*% result$par;
auto <- acf(ret);

Akaike <- matrix(nrow=3,ncol=4);
Beysian <- matrix(nrow=3,ncol=4);
for (p in 0:2) {
    for (q in 1:4) {
        model <- arima(ret, order=c(p,0,q));
        Akaike[p+1,q] <- model$aic;
        Beysian[p+1,q] <- -2*model$loglik + (p+q+3)*log(T);
    }
}
## MA(1) or ARMA(2,3)
n = T/4;
model <- Arima(ret[1:T-n], order=c(0,0,1));
predicted <- rep(NA, n);
for (i in 1 : n) {
    predicted[i] <- forecast.Arima(model, 1);
    if (i %% 10 == 0) {
        model <- Arima(ret[1:T-n], order=c(0,0,1));        
    }
}

