rm(list=ls());
graphics.off();

library("fGarch");
source("../../r/libxxie.r");

## currencies <- c(
##     "AUD_SEK_Rates",
##     "CAD_SEK_Rates",
##     ## "CHF_SEK_Rates",

##     "CNY_SEK_Rates",
##     "CZK_SEK_Rates",
##     "DKK_SEK_Rates",

##     "EUR_SEK_Rates",
##     "GBP_SEK_Rates",
##     "HKD_SEK_Rates",

##     "HUF_SEK_Rates",
##     "JPY_SEK_Rates",
##     "KRW_SEK_Rates",

##     "MAD_SEK_Rates",
##     "MXN_SEK_Rates",
##     "NOK_SEK_Rates",

##     "NZD_SEK_Rates",
##     ## "PLN_SEK_Rates",
##     ## "SAR_SEK_Rates",

##     "SGD_SEK_Rates",
##     ## "THB_SEK_Rates",
##     ## "TRY_SEK_Rates",

##     "USD_SEK_Rates"
##     );

## names <- c(
##     "AUD",
##     "CAD",
##     ## "CHF",

##     "CNY",
##     "CZK",
##     "DKK",

##     "EUR",
##     "GBP",
##     "HKD",

##     "HUF",
##     "JPY",
##     "KRW",

##     "MAD",
##     "MXN",
##     "NOK",

##     "NZD",
##     ##"SAR",
##     "SGD",
##     "USD"
##     );

assets = c("DAX", "FTSE100", "SP500", "DJIA", "OMXS30");
R <- getAssetReturns("2012-01-01", "2015-01-01", assets, 1, "closing", "localhost");

coef <- matrix(NA, nrow=length(assets), ncol=3);
for (i in 1:length(assets)) {
    X <- getAssetReturns("2012-01-01", "2015-01-01", assets[i], 1, "closing", "localhost");
    M <- garchFit(~garch(2, 1), data=X, trace=FALSE);
    coef[i, ] <- M@fit$params$params[c(3,4,7)];    
}

X <- getAssetReturns("2012-01-01", "2015-01-01", "SP500", 1, "closing", "localhost");
M <- garchFit(~garch(2, 1), data=X, trace=FALSE);


##            [,1]       [,2]      [,3]
## [1,] 0.02749864 0.04228535 0.8968533
## [2,] 0.10623464 0.02904907 0.7829784
## [3,] 0.08835834 0.09685783 0.6543018
## [4,] 0.07490443 0.12812118 0.6543123
## [5,] 0.01677130 0.03622593 0.9232396



## coef <- matrix(NA, nrow=p, ncol=3);
## for (i in 1:p) {
##     M <- garchFit(~garch(2,1), data=X[, i], trace=FALSE);
##     coef[i, ] <- M@fit$params$params[c(3,4,7)];
## }

