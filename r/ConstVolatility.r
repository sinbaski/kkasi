rm(list=ls());
graphics.off();
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

C <- cov(X);

Y <- mvrnorm(n=500*n, mu=rep(0, p), Sigma=C);

X2 <- X^2;
Y2 <- Y^2;

QX <- matrix(NA, p, p);
for (i in 1:p) {
    for (j in 1:i) {
        r <- quantile(X2[, i] * X2[, j], 1 - 1/n);
        QX[i, j] = r;
        QX[j, i] = r;
    }
}
CX <- cov(X2 - matrix(rep(apply(X2, MARGIN=2, FUN=mean), n), byrow=TRUE, n, p)) * dim(X2)[1] / max(QX);
E <- eigen(CX);

QY <- matrix(NA, p, p);
for (i in 1:p) {
    for (j in 1:i) {
        r <- quantile(Y2[, i] * Y2[, j], 1 - 1/dim(Y2)[1]);
        QY[i, j] = r;
        QY[j, i] = r;
    }
}
CY <- cov(Y2 - matrix(rep(apply(Y2, MARGIN=2, FUN=mean), n), byrow=TRUE, dim(Y2)[1], p)) * dim(Y2)[1] / max(QY);
D <- eigen(CY);

pdf("/tmp/iid_normal_eigenvalues.pdf");
plot(1:p, (D$values)/sum(D$values),
     main="FX and iid spectrum",
     ylim=c(0, 0.8),
     xlab=expression(i), ylab="", cex=2,
     pch=16, col="#FF0000");
points(1:p, (E$values)/sum(E$values), col="#000000", cex=2, pch=0);

legend("topright",
       legend=c(expression(cov(X^2)), expression(cov(Y^2))),
       col=c("#000000", "#FF0000"),
       pch=c(0, 16), cex=2);
grid();
dev.off();

pdf("/tmp/iid_normal_eigenvectors.pdf", width=20, height=10);
par(mfrow=c(3,6));
for (i in 1:p) {
    V <- E$vectors[, i];
    U <- D$vectors[, i];
    
    if (sum(abs(V - U)) > sum(abs(V + U))) {
        U <- -U;
    }
    s <- sign(V[which.max(abs(V))]);

    ## k <- which.max(abs(V));
    ## V <- V * sign(V[k]);
    plot(1:p, V * s, main=sprintf("FX & iid V[%d]", i),
         xlab="i", ylab=expression(V[i]),
         ylim=c(-1, 1), pch=0,
         xaxt="n");
    axis(side=1, at=1:p, labels=names, las=2);
    
    ## if (!(i %in% c(2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16))) {
    ##     V <- D$vectors[, i];
    ##     k <- which.max(abs(V));
    ##     V <- V * sign(V[k]);
    ## }
    points(1:p, U * s, main=sprintf("eigenvector[%d]", i),
           col="#FF0000", pch=16);

    ## V <- F$vectors[, i];
    ## k <- which.max(abs(V));
    ## V <- V * sign(V[k]);
    ## points(1:p, V, main=sprintf("eigenvector[%d]", i),
    ##        col="#00FF00", pch=17);
    grid();
}
dev.off();


