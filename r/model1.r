library(RMySQL);
rm(list=ls());

n_obs = 1000;
database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host='localhost');
R = matrix(nrow=n_obs, ncol=3);
ptfl = data.frame(c("ATCO_series_A_ST_SE",
    "ATCO_series_B_ST_SE",
    "INVE_series_B_ST_SE"),
    c(0.7568, -0.6412, -0.0847));
for (k in 1:nrow(ptfl)) {
    results = dbSendQuery(database,
        sprintf("select closing from
         %s order by day desc limit %d;", ptfl[k,1], n_obs+1));
    prices = rev(log(fetch(results, n=-1)[[1]]));
    dbClearResult(results);
    
    R[,k] = tail(prices, -1) - head(prices, -1);
}

dbDisconnect(database);

X = R %*% ptfl[,2];
Y = X;
Akaike <- matrix(nrow=5,ncol=5);
Beysian <- matrix(nrow=5,ncol=5);
for (p in 1:5) {
    for (q in 1:5) {
        model = arima(Y, order=c(p,0,q));
        Akaike[p,q] <- model$aic;
        Beysian[p,q] <- -2*model$loglik + (p+q+1)*log(length(Y));
    }
}
    
## Y <- diff(X, differences=1, lags=3);
models = arima(Y, order=c(2, 0, 2));
acf(models$residuals);
