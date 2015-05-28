library(RMySQL);
rm(list=ls());

source("libxxie.r");

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
models <- matrix(NA, nrow=p+1, ncol=p);
par(mfrow=c(3,4));
for (k in 1:dim(M)[2]) {
    mdl <- lm(R[,k] ~ factors);
    models[,k] <- mdl$coef / sum(abs(mdl$coef));
    acf(mdl$residuals, main=tables[k]);
}

A = R[, 10] - factors %*% models[-1,10];



    
