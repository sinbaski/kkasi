library(RMySQL);
rm(list=ls());

n_obs = 1000;
database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host=Sys.getenv("PB"));
## ptfl = data.frame(c("ELUX_series_B_ST_SE",
##     "ERIC_series_B_ST_SE",
##     "LUPE_ST_SE",
##     "NDA_series_SEK_ST_SE",
##     "SKA_series_B_ST_SE"),
##     c(0.4306, 0.0656, -0.1065, 0.2999, 0.0974));
ptfl = data.frame(c("ATCO_series_A_ST_SE",
    "ATCO_series_B_ST_SE"),
    c(0.7568, -0.7000));
R = matrix(nrow=n_obs, ncol=dim(ptfl)[1]);
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
    

bestModel = arima(Y, order=c(4,0,4));
