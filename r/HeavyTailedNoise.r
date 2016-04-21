rm(list=ls());
## library(fields);
source("libxxie.r");

#Simulate samples described by a heavy-tailed SV model.
## n <- 100000;
## p <- 20;
## X <- matrix(rt(n=n*p, df=3), nrow=n, ncol=p);
## Y <- matrix(NA, n, p);
## eta <- matrix(rnorm(n=(p+1)*n), nrow=n, ncol=p+1);
## for (i in 2:(p+1)) {
##     Y[, i-1] <- X[, i-1] * exp(eta[, i-1] + eta[,i]);
## }

## Use FX data instead.
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
Y <- getAssetReturns("2010-01-04", "2016-04-01", currencies, 1, "rate", "31.208.142.23");
p <- dim(Y)[2];
E <- eigen(cov(Y));

pdf("~/Documents/FX_eigenvalues.pdf")
plot(1:p, E$values, xlab="i", ylab=expression(lambda[(i)]));
dev.off();

pdf("~/Documents/FX_ratios.pdf");
plot(1:(p-1), E$values[2:p]/E$values[1:(p-1)],
     type="b", ylim=c(0, 1),
     xlab="i", ylab=expression(lambda[(i+1)]/lambda[(i)]));
dev.off();

pdf("~/Documents/FX_eigenvector1.pdf")
plot(1:p, E$vectors[, 1], xlab="i", ylab=expression(V["i,1"]));
dev.off();


## tail.indices <- matrix(NA, p, p);
## for (i in (1:p)) {
##     for (j in (i:p)) {
##         tail.indices[i,j] <- pikandsEstimate(Y[,i] * Y[,j], prob=0.97);
##         tail.indices[j,i] <- tail.indices[i,j]
##     }
## }

## M = max(tail.indices);
## m <- min(tail.indices);

## colors <- gray((tail.indices - m)/(M - m));
## pdf("/tmp/ParetoNoise_AbsnormVolatility_Hill.pdf")
## plot(1, 1, type="n", xlim=c(1,p+2), ylim=c(1,p));
## for (i in 1:p) {
##     for (j in i:p) {
##         points(x=i, y=j, pch=19, cex=3, col=colors[(j-1)*p+i]);
##         ## points(x=i, y=j, pch=19, cex=3,
##         ##        col=colors[floor((tail.indices[i,j]-m)/(M-m)*length(colors))+1]);
##     }
## }
## ## axis(side=1, at=1:p, labels=names);
## ## axis(side=2, at=1:p, labels=names);
## ## image.plot(legend.only=TRUE, zlim=c(m, M), col=sort(colors))
## image.plot(x=seq(1,p+2), y=1:p, legend.only=TRUE, zlim=c(m, M), col=sort(colors));
## dev.off();

