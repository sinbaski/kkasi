rm(list=ls());
library(abind);
library(RMySQL);
library(fields);
source("libxxie.r");

currencies <- c(
    ## Oceania
    "NZD_SEK_Rates",
    "AUD_SEK_Rates",
    
    ## Asia
    "CNY_SEK_Rates",
    "HKD_SEK_Rates",
    "JPY_SEK_Rates",
    "KRW_SEK_Rates",
#    "SAR_SEK_Rates", # Saudi-Arabia
    "SGD_SEK_Rates", # Singapore
    "THB_SEK_Rates", # Thailand
#    "TRY_SEK_Rates", # Turkey
    
    ## Europe
    "CHF_SEK_Rates",
#    "CZK_SEK_Rates", # Czech
    "DKK_SEK_Rates",
    "EUR_SEK_Rates",
    "GBP_SEK_Rates",
#    "HUF_SEK_Rates", # Hungary
    "NOK_SEK_Rates",
#    "PLN_SEK_Rates", # Poland

    ## Africa
#    "MAD_SEK_Rates", # Maroco
    
    ## Americas
    "CAD_SEK_Rates",
    "USD_SEK_Rates",
    "MXN_SEK_Rates" # Mexico
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
            "select rate from %s order by day;",
            currencies[i])
        );
    prices <- fetch(results, n=-1)[[1]];
    R <- diff(log(prices), lag=1);
    ## R <- R - mean(R);
    if (length(ret) == 0) {
        ret <- R;
    } else {
        ret <- cbind(ret, R);
    }
    dbClearResult(results);
}
dbDisconnect(database);

test <- rep(NA, 21);
m <- 1;
for (n in seq(from=200, to=1200, by=50)) {
    C <- cov(ret[1:n,]);
    E <- eigen(C, only.values=TRUE);
    test[m] <- (E$values[1] - E$values[p])/sum(E$values);
    m <- m + 1;
}
pdf("/tmp/SV.pdf");
plot(seq(from=200, to=1200, by=50), test, type="p",
     main=expression(frac(lambda[(1)] - lambda[(p)], plain(trace))),
     xlab="n", ylab="");
dev.off();
