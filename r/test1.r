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
auto <- acf(ret);

Akaike <- matrix(nrow=5,ncol=5);
Beysian <- matrix(nrow=5,ncol=5);
for (p in 1:5) {
    for (q in 1:5) {
        model <- Arima(ret, order=c(p,0,q));
        Akaike[p,q] <- model$aic;
        Beysian[p,q] <- model$bic;
    }
}

## ARMA(3,1)
n = floor((T-1)/4);
model <- Arima(ret[1:(T-n)], order=c(3,0,1));

predicted <- rep(NA, n);
for (i in 1 : n) {
    predicted[i] <- forecast.Arima(model, h=1)$mean;
    model <- Arima(ret[1:(T-n+i)], model=model);
}

