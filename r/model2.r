library(RMySQL);
library(alabama);
rm(list=ls());

source("libxxie.r");

my.fun <- function(par, Y, X) {
    return(sum((Y - X %*% par)^2));
}

my.fun.der <- function(par, Y, X) {
    return(-2*t(X) %*% Y + 2*t(X) %*% X %*% par);
}

heq <- function(par, Y, X) {
    return(sum(abs(par)) - 1);
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

C <- t(R) %*% R / dim(R)[1];
E <- eigen(C, symmetric=TRUE);
M <- E$vectors;
for (k in 1 : dim(M)[2]) {
    M[,k] <- M[,k]/sum(abs(M[,k]));
}
factors <- R %*% M;

par(mfrow=c(3,4));
for (k in 1:p) {
    A <- rep(0, p);
    A[k] <- 1;
    result <- auglag(par=A, fn=my.fun, gr=my.fun.der, Y=R[,k], X=factors,
                     heq=heq);
    res = R[,k] - factors %*% result$par;
    acf(res, main=tables[k]);
}

ind <- rep(0,p);
for (k in 1 : p) {
    if (max(abs(M[,k] - abs(M[,k]))) < 1.0e-4) {
        ind[k] <- 1;
    }
}
