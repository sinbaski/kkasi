rm(list=ls());
graphics.off();
library("tseries");
library("fGarch");
library("MASS");
source("libxxie.r");

currencies <- c(
    "AUD_SEK_Rates",
    "CAD_SEK_Rates",
    ## "CHF_SEK_Rates",

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
    ## "PLN_SEK_Rates",
    ## "SAR_SEK_Rates",

    "SGD_SEK_Rates",
    ## "THB_SEK_Rates",
    ## "TRY_SEK_Rates",

    "USD_SEK_Rates"
    );

names <- c(
    "AUD",
    "CAD",
    ## "CHF",

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
    ##"SAR",
    "SGD",
    "USD"
    );

X <- getAssetReturns("2010-01-04", "2016-04-01", currencies, 1,
                     "rate", "localhost");
n <- dim(X)[1];
p <- dim(X)[2];

episode <- 20;
lam1 <- rep(NA, floor(n/episode));
for (k in 1:length(lam1)) {
    Y <- X[((k-1)*episode + 1) : (k * episode), ];
    ## for (i in 1:p) {
    ##     M <- garchFit(~garch(1,1), data=X[, i], trace=FALSE);
    ##     coef[i, ] <- M@fit$params$params[c(2,3,5)];
    ##     res[, i] <- M@residuals;
    ##     C <- cor(res);
    ## }
    CY <- cov(Y);
    E <- eigen(CY);
    lam1[k] <- E$values[1];
}
pdf("/tmp/Lambda1_evolution.pdf")
plot(1:length(lam1), lam1, type="l",
     xlab=expression(k), ylab=expression(lambda[1]^{(1)}));
dev.off();

