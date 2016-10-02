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

inno <- matrix(NA, nrow=n, ncol=p);
coef <- matrix(NA, nrow=p, ncol=3);
ics <- matrix(NA, nrow=p, ncol=4);
for (i in 1:p) {
    ## spec <- ugarchspec(
    ##     mean.model=list(
    ##         armaOrder=c(0,0),
    ##         include.mean=FALSE
    ##     ),
    ##     variance.model=list(
    ##         garchOrder=c(1, 1)
    ##     ),
    ##     distribution.model="norm"
    ## );

    ## solvers <- c("nlminb", "solnp", "lbfgs", "gosolnp", "nloptr", "hybrid");
    ## M1 <- NULL;
    ## j <- 1;
    ## while (j <= length(solvers) && (is.null(M1) || convergence(M1) != 0)) {
    ##     M1 <- ugarchfit(data=X[, i] - mean(X[, i]), spec=spec, solver=solvers[j]);
    ##     if (convergence(M1) != 0) {
    ##         j <- j + 1;
    ##         next
    ##     }
    ##     inno[, i] <- residuals(M1) / sigma(M1);
    ##     inno[, i] <- inno[, i] - mean(inno[, i]);
    ##     coef[i, ] <- coef(M1);
    ## }
    ## if (is.null(M1) || convergence(M1)) {
    ##     stop(sprintf("Estimation failed for %d", i));
    ## }
    
    M2 <- garchFit(~garch(1,0),
                  data=X[, i],
                  trace=FALSE,
                  cond.dist="norm",
                  ## shape=4,
                  include.shape=FALSE,
                  include.mean=TRUE,
                  include.delta=FALSE,
                  include.skew=FALSE
                  );
    ## coef[i, ] <- M2@fit$params$params[c(2,3,5)];
    coef[i, ] <- c(M2@fit$coef[-1], 0);
    ## coef[i, ] <- M2@fit$coef[-1];
    inno[, i] <- M2@residuals / M2@sigma.t;
    inno[, i] <- inno[, i] - mean(inno[, i]);
    inno[, i] <- inno[, i] / sd(inno[, i]);
    ics[i, ] <- M2@fit$ics;
}
write.table(x=ics, file="GARCH_ic.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);
## mean.inno <- apply(inno, MARGIN=2, FUN=mean);
## inno <- inno - matrix(rep(mean.inno, n), nrow=n, ncol=p, byrow=TRUE);
## inno <- inno %*% diag(1 / apply(inno, MARGIN=2, FUN=sd));
C <- cor(inno);
## V <- apply(res, MARGIN=2, FUN=sd);
## res <- res / matrix(rep(V, n), byrow=TRUE, nrow=n, ncol=p);
## sigma <- X / res;
## M <- t(sigma) %*% sigma / n;

## U <- eigen(M * C);
## vol <- matrix(NA, nrow=n, ncol=p);
## for (i in 1:p) {
##     vol[1, i] <- coef[i, 1];
##     for (t in 2:(n-10)) {
##         vol[t, i] <- coef[i, 1] + coef[i, 1] * X[t-1, i]^2 + coef[i, 2] * vol[t-1, i];
##     }
## }
## vol <- sqrt(vol);

## V <- t(vol) %*% vol / n;
## C <- cov(res);

## diag(C) <- 1;
## write.table(x=cor(res), file="/tmp/cor.txt", row.names=FALSE, col.names=FALSE);
W <- matrix(NA, nrow=100*n, ncol=p);
sig2 <- matrix(NA, nrow=dim(W)[1], ncol=p);
# set the initial values
for (i in 1:p) {
    sig2[1, i] <- coef[i, 1]/(1 - coef[i, 3]);
    ## sig2[1, i] <- 0;
}
for (i in 1:dim(W)[1]) {
    eta <- rmvnorm(n=1, mean=rep(0, p), sigma=C);
    ## eta <- rmvt(n=1, sigma=C, df=4);
    W[i, ] <- eta * sqrt(sig2[i,]);
    if (i < dim(W)[1])
        sig2[i+1, ] <- coef[, 2] * W[i, ]^2 + coef[, 3] * sig2[i, ] + coef[, 1];
}
## U <- sqrt(sig2);
## S <- t(U) %*% U / dim(sig2)[1];
## Compute Hill Estimators
# tailIndices <- matrix(NA, p, p);
## tailIndices <- rep(NA, p);
## for (i in 1:p) {
##     a <- hillEstimate(X[, i]^2, prob=0.97);
##     tailIndices[i] <- a;
##     ## for (j in 1:i) {
##     ##     T <- W[,i] * W[,j];
##     ##     a <- hillEstimate(T, prob=0.97);
##     ##     # tailIndices[(j-1)*j/2 + i] <- a;
##     ##     tailIndices[i, j] <- a;
##     ## }
## }

## write.table(format(tailIndices, digits=2),
##             quote=FALSE, sep="  ",
##             row.names=FALSE, col.names=FALSE,
##             file="/tmp/Hill_Simulated.txt");

## pdf("/tmp/FX_real_n_simulated_eigenvalues.pdf", width=14, height=14);
## M <- apply(X, MARGIN=2, FUN=mean);
## data <- X - matrix(rep(M, n), nrow=n, ncol=p, byrow=TRUE);
# q <- max(apply(data^2, MARGIN=2, FUN=quantile, probs=0.99));

## QX <- matrix(NA, p, p);
## for (i in 1:p) {
##     for (j in 1:i) {
##         r <- quantile(X2[, i] * X2[, j], 1 - 1/n);
##         QX[i, j] = r;
##         QX[j, i] = r;
##     }
## }

## CX <- cov(X2 - matrix(rep(apply(X2, MARGIN=2, FUN=mean), n), byrow=TRUE, n, p)) * dim(X2)[1] / max(QX);

## CX <- cov(X2) * dim(X2)[1] / max(QX);
CX <- cov(X);
E <- eigen(CX);

CY <- cov(inno);
D <- eigen(CY);

## QW <- matrix(NA, p, p);
## for (i in 1:p) {
##     for (j in 1:i) {
##         r <- quantile(W2[, i] * W2[, j], 1 - 1/dim(W2)[1]);
##         QW[i, j] = r;
##         QW[j, i] = r;
##     }
## }
## CW <- cov(W2 - matrix(rep(apply(W2, MARGIN=2, FUN=mean), n), byrow=TRUE, dim(W2)[1], p)) * dim(W2)[1] / max(QW);
# CW <- cov(W2) * dim(W2)[1] / max(QW);
CW <- cov(W);
F <- eigen(CW);

pdf("/tmp/FX_eigenvalues.pdf");
## plot(1:p, sig.eig$values, type="p", pch=17,
##      main="FX and GARCH(1,1) spectrum", col="#00FF00"
## );
plot(1:p, E$values/sum(E$values), type="p", pch=0,
     main="FX and ARCH(1) spectrum"
);
points(1:p, (D$values)/sum(D$values), col="#0000FF", pch=17);
points(1:p, (F$values)/sum(F$values), pch=16, col="#FF0000");

## ## points(1:p, (E1$values)/sum(E1$values), col="#FF0000", cex=2, pch=15);
## ## points(1:p, (F1$values)/sum(F1$values), col="#00FF00", cex=2, pch=16);
## ## points(1:p, (F1$values)/sum(F1$values), col="#00FF00", cex=2, pch=17);

legend("topright",
##       legend=c(expression(sigma[i] * sigma[j]), expression(cov(W)), expression(cov(X))),
       legend=c(expression(cov(FX)), expression(cov(inno)), expression(cov(sim.))),
       col=c("#000000", "#0000FF", "#FF0000"),
       pch=c(0, 17, 16));
grid();
dev.off();

## pdf("/tmp/first_eigenvector.pdf");
## V <- sig.eig$vectors[, 1];
## V <- V * sign(V[which.max(abs(V))])
## plot(1:p, V, type="p", xaxt="n", xlab="i", ylab=expression(V[1](i)));
## axis(side=1, at=1:p, labels=names, las=2);
## dev.off();


pdf("/tmp/FX_eigenvectors.pdf", width=20, height=10);
par(mfrow=c(3,6));
mse <- c(0, 0);
for (i in 1:p) {
    V <- E$vectors[, i];
    U <- D$vectors[, i];
    Q <- F$vectors[, i];
    
    if (sum(abs(V - U)) > sum(abs(V + U))) {
        U <- -U;
    }
    if (sum(abs(V - Q)) > sum(abs(V + Q))) {
        Q <- -Q;
    }
    s <- sign(V[which.max(abs(V))]);

    mse[1] <- mse[1] + sum(abs(V * s - U * s));
    mse[2] <- mse[2] + sum(abs(V * s - Q * s));
    
    plot(1:p, V * s, main=sprintf("FX & ARCH(1) V[%d]", i),
         xlab="i", ylab=expression(V[i]),
         ylim=c(-1, 1), pch=0,
         xaxt="n");
    axis(side=1, at=1:p, labels=names, las=2);
    
    ## points(1:p, U * s,
    ##        col="#0000FF", pch=17);
    points(1:p, Q * s, col="#FF0000", pch=16);

    grid();
}
dev.off();
