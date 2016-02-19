rm(list=ls());
library(abind);
library(RMySQL);
source("libxxie.r");

currencies <- c(
    "AUD_SEK_Rates",
    "CAD_SEK_Rates",
    "CHF_SEK_Rates",

    "CNY_SEK_Rates",
    "CZK_SEK_Rates",
    "DKK_SEK_Rates",

    "EUR_SEK_Rates",
    "GBP_SEK_Rates",
    "HKD_SEK_Rates",

    "HUF_SEK_Rates",
    "JPY_SEK_Rates",
    "KRW_SEK_Rates",

    "MAD_SEK_Rates",
    "MXN_SEK_Rates",
    "NOK_SEK_Rates",

    "NZD_SEK_Rates",
    "PLN_SEK_Rates",
    "SAR_SEK_Rates",

    "SGD_SEK_Rates",
    "THB_SEK_Rates",
    "TRY_SEK_Rates",

    "USD_SEK_Rates"
    );


p <- length(currencies);
tailIndices <- rep(0, choose(p,2));
ret <- list();

database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host="localhost");
for (i in 1:p) {
    results <- dbSendQuery(
        database,
        sprintf(
            "select rate from %s where day >= '2013-01-01' order by day;",
            currencies[i])
        );
    prices <- fetch(results, n=-1)[[1]];
    R <- diff(log(prices), lag=1);
    R <- R - mean(R);
    if (length(ret) == 0) {
        ret <- R;
    } else {
        ret <- cbind(ret, R);
    }
    dbClearResult(results);
}
dbDisconnect(database);

for (i in 1:p) {
    for (j in i:p) {
        X <- ret[,i] * ret[,j];
        tailIndices[(j-1)*j/2 + i] <- hillEstimate(X, probs=0.97);
##        tailIndices[j,i] <- tailIndices[i,j];
    }
}

## filled.contour(x=1:p, y=1:p, z=tailIndices, axes = TRUE);
names <- c(
    "AUD",
    "CAD",
    "CHF",
    "CNY",
    "CZK",
    "DKK",
    "EUR",
    "GBP",
    "HKD",
    "HUF",
    "JPY",
    "KRW",
    "MAD",
    "MXN",
    "NOK",
    "NZD",
    "PLN",
    "SAR",
    "SGD",
    "THB",
    "TRY",
    "USD"
    );
M <- max(tailIndices);
m <- min(tailIndices);
colors <- gray((tailIndices-m)/(M-m));
pdf("/tmp/FX_HillEstimates.pdf")
plot(1, 1, type="n", xlim=c(1,p), ylim=c(1,p));
for (i in 1:p) {
    for (j in i:p) {
        points(x=i, y=j, pch=19, cex=(M/tailIndices[j*(j-1)/2 + i])^2,
               col=colors[(j-1)*j/2+i]);
    }
}
axis(side=1, at=1:p, labels=names);
axis(side=2, at=1:p, labels=names);
dev.off();

