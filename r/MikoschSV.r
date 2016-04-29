rm(list=ls());
library(abind);
library(RMySQL);
library(fields);
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
    ## "PLN_SEK_Rates",
    ## "SAR_SEK_Rates",

    "SGD_SEK_Rates",
    ## "THB_SEK_Rates",
    ## "TRY_SEK_Rates",

    "USD_SEK_Rates"
    );

ret <- getAssetReturns("2010-01-04", "2016-04-01",
                       currencies, 1,
                       "rate", "localhost");
n <- dim(ret)[1];
p <- length(currencies);
tailIndices <- rep(0, choose(p,2));

for (i in 1:p) {
    for (j in i:p) {
        X <- ret[,i] * ret[,j];
        tailIndices[(j-1)*j/2 + i] <- hillEstimate(X, prob=0.97);
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
    "SGD",
    "USD"
    );

M <- max(tailIndices);
m <- min(tailIndices);
colors <- gray((tailIndices-m)/(M-m));
pdf("/tmp/FX_HillEstimates.pdf")
plot(1:p, 1:p, type="n", xlim=c(1,p+2), ylim=c(1,p),
     xaxt='n', yaxt='n', xlab="", ylab="");
for (i in 1:p) {
    for (j in i:p) {
        ## points(x=i, y=j, pch=19, cex=(M/tailIndices[j*(j-1)/2 + i])^2,
        ##        col=colors[(j-1)*j/2+i]);
        points(x=i, y=j, pch=19, cex=3.5, col=colors[(j-1)*j/2+i]);
    }
}
## grid();
axis(side=1, at=1:p, labels=names, las=2);
axis(side=2, at=1:p, labels=names, las=2);
## image.plot(legend.only=TRUE, zlim=c(m, M), col=sort(colors))
image.plot(x=seq(1,p+2), y=1:p, legend.only=TRUE, zlim=c(m, M), col=sort(colors));
dev.off();

