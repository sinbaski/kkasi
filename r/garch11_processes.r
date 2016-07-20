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

# res <- matrix(NA, nrow=n-1, ncol=p);
res <- matrix(NA, nrow=n, ncol=p);
coef <- matrix(NA, nrow=p, ncol=3);
for (i in 1:p) {
    ## M <- garch(x=X[, i], order=c(1, 1), trace=FALSE);
    ## coef[i, ] <- M$coef;
    ## res[, i] <- tail(M$residuals, -1);

    M <- garchFit(~garch(1,1), data=X[, i], trace=FALSE);
    coef[i, ] <- M@fit$params$params[c(2,3,5)];
    res[, i] <- M@residuals;
    
    ## M <- estimGARCH(0, 0.01, 0, X[, i]);
    ## coef[i, ] <- M$coef;
    ## res[, i] <- M$residus;
    ## print(c(names[i], coef[i, 2:3]));
}

## C <- cov(res);
C <- cor(res);
## diag(C) <- 1;
write.table(x=cor(res), file="/tmp/cor.txt", row.names=FALSE, col.names=FALSE);
Y <- matrix(NA, nrow=100*n, ncol=p);
sig2 <- matrix(NA, nrow=dim(Y)[1], ncol=p);
# set the initial values
for (i in 1:p) {
    sig2[1, i] <- coef[i, 1]/(1 - coef[i, 3]);
}
for (i in 1:dim(Y)[1]) {
    eta <- mvrnorm(n=1, mu=rep(0, p), Sigma=C);
    Y[i, ] <- eta * sqrt(sig2[i,]);
    if (i < dim(Y)[1])
        sig2[i+1, ] <- coef[, 2] * Y[i, ]^2 + coef[, 3] * sig2[i, ] + coef[, 1];
}

## Compute Hill Estimators
tailIndices <- matrix(NA, p, p);
for (i in 1:p) {
    for (j in 1:i) {
        T <- Y[,i] * Y[,j];
        a <- hillEstimate(T, prob=0.95);
        # tailIndices[(j-1)*j/2 + i] <- a;
        tailIndices[i, j] <- a;
    }
}

write.table(format(tailIndices, digits=2),
            quote=FALSE, sep="  ",
            row.names=FALSE, col.names=FALSE,
            file="/tmp/Hill_Simulated.txt");

## pdf("/tmp/FX_real_n_simulated_eigenvalues.pdf", width=14, height=14);
M <- apply(X, MARGIN=2, FUN=mean);
Q <- apply(X, MARGIN=2, quantile, 1-1/n);
data <- X - matrix(rep(M, n), nrow=n, ncol=p, byrow=TRUE);
C <- (t(data) %*% data)/max(Q)^2;
#E <- eigen(cov(X - mean(X)));
E <- eigen(C);

M <- apply(Y, MARGIN=2, FUN=mean);
Q <- apply(Y, MARGIN=2, quantile, 1-1/dim(Y)[1]);
data <- Y - matrix(rep(M, dim(Y)[1]), nrow=dim(Y)[1], ncol=p, byrow=TRUE);
C <- (t(data) %*% data)/max(Q)^2;
D <- eigen(C);

plot(1:p, (E$values)/sum(E$values),
     main=expression(lambda[(i)]/trace),
     ylim=c(0, 0.8),
     xlab=expression(i), ylab="", cex=2, pch=0);
points(1:p, (D$values)/sum(D$values), col="#FF0000", cex=2, pch=16);

## ## points(1:p, (E1$values)/sum(E1$values), col="#FF0000", cex=2, pch=15);
## ## points(1:p, (D1$values)/sum(D1$values), col="#00FF00", cex=2, pch=16);
## ## points(1:p, (F1$values)/sum(F1$values), col="#00FF00", cex=2, pch=17);

legend("topright",
       legend=c(expression(X - EX), expression(Y - EY)),
       col=c("#000000", "#FF0000"),
       pch=c(0, 16), cex=2);
grid();
## dev.off();

## pdf("/tmp/FX_real_n_simulated_eigenvectors.pdf", width=20, height=10);
par(mfrow=c(3,6));
for (i in 1:p) {
    V <- E$vectors[, i];
    k <- which.max(abs(V));
    V <- V * sign(V[k]);
    plot(1:p, V, main=sprintf("eigenvector[%d]", i),
         xlab="i", ylab=expression(V[i]),
         ylim=c(-1, 1), pch=0,
         xaxt="n");
    axis(side=1, at=1:p, labels=names, las=2);
    
    if (!(i %in% c(2, 6, 8))) {
        V <- D$vectors[, i];
        k <- which.max(abs(V));
        V <- V * sign(V[k]);
    }
    points(1:p, V, main=sprintf("eigenvector[%d]", i),
           col="#FF0000", pch=16);

    ## V <- F$vectors[, i];
    ## k <- which.max(abs(V));
    ## V <- V * sign(V[k]);
    ## points(1:p, V, main=sprintf("eigenvector[%d]", i),
    ##        col="#00FF00", pch=17);
    grid();
}
## dev.off();
