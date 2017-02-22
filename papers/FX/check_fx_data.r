graphics.off();
database = dbConnect(MySQL(), user='root', password='q1w2e3r4',
                     dbname='avanza', host="localhost");

result <- dbSendQuery(database,
                      "select distinct day from USD_SEK_Rates where day between \"2010-01-04\" and \"2016-04-01\";")
D = fetch(result, n=-1)[[1]];
X <- getAssetPrices("2010-01-04", "2016-04-01", currencies, 1,
                     "rate", "localhost");

par(mfrow=c(3,6));
## mse <- c(0, 0);
for (i in 1:p) {
    plot(as.xts(X[, i], order.by=as.Date(D)),
         main=names[i]);
}

