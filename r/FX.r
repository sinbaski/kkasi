rm(list=ls());
graphics.off();
## require("tseries");
require("fGarch");
require("rugarch");
require("mvtnorm");
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
X <- X - matrix(rep(apply(X, MARGIN=2, FUN=mean), n), nrow=n, ncol=p, byrow=TRUE);
N <- p;

E <- cov(X);
Eig <- eigen(E);
Y <- X %*% Eig$vectors[, 1:N];
D <- cov(Y);
Dig <- eigen(D);
    
C <- matrix(NA, nrow=N, ncol=p);
Res <- matrix(NA, nrow=n, ncol=p);
for (i in 1:p) {
    R <- X[, i];
    model <- lm(R ~ Y[, 1:N]);
    C[, i] <- coef(model)[1:N+1];
    Res[, i] <- residuals(model);
}
Dep <- cov(Res);
F <- t(C) %*% D %*% C;
Fig <- eigen(F);
