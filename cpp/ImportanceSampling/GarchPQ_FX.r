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

## assets = c("DAX", "FTSE100", "SP500", "DJIA", "OMXS30");
## R <- getAssetReturns("2012-01-01", "2015-01-01", assets, 1, "closing", "localhost");

## coef <- matrix(NA, nrow=length(assets), ncol=3);
## for (i in 1:length(assets)) {
##     X <- getAssetReturns("2012-01-01", "2015-01-01", assets[i], 1, "closing", "localhost");
##     garch21 <- garchFit(~garch(2, 1), data=X, trace=FALSE,
##                   include.mean=FALSE);
##     garch11 <- garchFit(~garch(1, 1), data=X, trace=FALSE,
##                   include.mean=FALSE);
##     coef[i, ] <- garch21@fit$coef[-1];
## }

X <- getAssetReturns("2012-01-01", "2015-01-01", "DJIA", 1, "closing", "localhost");
M21 <- garchFit(~garch(2, 1), data=X, trace=FALSE);
M11 <- garchFit(~garch(1, 1), data=X, trace=FALSE);
sum(M21@fit$coef[c(3,4,5)])
ics = matrix(c(M11@fit$ics, M21@fit$ics), byrow=T, nrow=2, ncol=4)


##           alpha1       alpha2     beta1   sum
## DAX     0.02074733 0.04104947 0.9102376   0.97
## FTSE100 0.08806776 0.04548650 0.7808275   0.91
## SP500   0.07949678 0.08765884 0.6683833   0.84
## DJIA    0.06157793 0.12795424 0.6610499   0.85
## OMXS30  0.00924104 0.04309249 0.9246996   0.98



## coef <- matrix(NA, nrow=p, ncol=3);
## for (i in 1:p) {
##     M <- garchFit(~garch(2,1), data=X[, i], trace=FALSE);
##     coef[i, ] <- M@fit$params$params[c(3,4,7)];
## }

