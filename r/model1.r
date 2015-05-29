library(RMySQL);
rm(list=ls());

n.obs = 1000;
database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host=Sys.getenv("PB"));
I = c(4, 5, 12)
ptfl = data.frame(c("ATCO_series_A_ST_SE",
    "ATCO_series_B_ST_SE",
    "INVE_series_B_ST_SE"),
    c(1, -0.8472, -0.1119));

results = dbSendQuery(database, sprintf("select symbol from %s;", 'OMXS30_components;'));
tables <- fetch(results, n=-1)[[1]];
tables <- gsub("[.]", "_", tables);
tables <- gsub("-", "_series_", tables);
tables <- sprintf("%s_SE", tables);
dbClearResult(results);

p = length(tables);

R = matrix(nrow=n.obs, ncol=length(tables));
for (k in 1:dim(R)[2]) {
    results = dbSendQuery(database,
        sprintf("select closing from
         %s order by day desc limit %d;", tables[k], n.obs+1));
    prices = rev(log(fetch(results, n=-1)[[1]]));
    dbClearResult(results);
    R[,k] = tail(prices, -1) - head(prices, -1);
}
dbDisconnect(database);

C <- t(R) %*% R / n.obs;
E <- eigen(C);
factors <- R %*% E$vectors;
par(mfrow=c(5,6));
for (k in 1:p) {
    ## mdl <- lm(R[,k] ~ factors);
    ## models[,k] <- mdl$coef[-1];
    ## acf(mdl$residuals, main=tables[k]);
    acf(factors[,k]);
}

model <- lm(R[,4] ~ factors[,1:20]);
mean(abs(model$residuals))
acf(model$residuals)

subject <- R[,I] %*% ptfl[[2]];


## X = R %*% ptfl[,2];
## Y = X;
## Akaike <- matrix(nrow=5,ncol=5);
## Beysian <- matrix(nrow=5,ncol=5);
## for (p in 1:5) {
##     for (q in 1:5) {
##         model = arima(Y, order=c(p,0,q));
##         Akaike[p,q] <- model$aic;
##         Beysian[p,q] <- -2*model$loglik + (p+q+1)*log(length(Y));
##     }
## }
    

## bestModel = arima(Y, order=c(4,0,4));
