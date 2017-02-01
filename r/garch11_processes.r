rm(list=ls());
graphics.off();
## require("tseries");
require("fGarch");
require("mvtnorm");
source("libxxie.r");

currencies <- c(
    "AUD_SEK_Rates",
    "CAD_SEK_Rates",
    ## "CHF_SEK_Rates",

    ## "CNY_SEK_Rates",
    "CZK_SEK_Rates",
    "DKK_SEK_Rates",

    "EUR_SEK_Rates",
    "GBP_SEK_Rates",
    "HKD_SEK_Rates",

    "HUF_SEK_Rates",
    "JPY_SEK_Rates",
    "KRW_SEK_Rates",

    ## "MAD_SEK_Rates",
    ## "MXN_SEK_Rates",
    "NOK_SEK_Rates",

    "NZD_SEK_Rates",
    "SGD_SEK_Rates",
    "USD_SEK_Rates"
    );

names <- c(
    "AUD",
    "CAD",
    ## "CHF",

    ## "CNY",
    "CZK",
    "DKK",

    "EUR",
    "GBP",
    "HKD",

    "HUF",
    "JPY",
    "KRW",

    ## "MAD",
    ## "MXN",
    "NOK",

    "NZD",
    "SGD",
    "USD"
    );

X <- getAssetReturns("2010-01-04", "2016-04-01", currencies, 1,
                     "rate", "localhost");
## X <- getInterpolatedReturns("2008-01-01", "2015-01-06",
## "DAX_components");

n <- dim(X)[1];
p <- dim(X)[2];
X <- X - matrix(rep(apply(X, FUN=mean, MARGIN=2), n), nrow=n, ncol=p, byrow=TRUE);

inno <- matrix(NA, nrow=n, ncol=p);
vola <- matrix(NA, nrow=n, ncol=p);
params <- matrix(0, nrow=p, ncol=3);
## ARCH
## params <- matrix(0, nrow=p, ncol=2);
ics <- matrix(NA, nrow=p, ncol=4);
for (i in 1:p) {
    M <- garchFit(~garch(1,1),
                  data=X[, i],
                  trace=FALSE,
                  cond.dist="norm",
                  ## shape=4,
                  include.shape=FALSE,
                  include.mean=FALSE,
                  include.delta=FALSE,
                  include.skew=FALSE
                  );
    params[i, ] <- coef(M);
    inno[, i] <- M@residuals / M@sigma.t;
    inno[, i] <- inno[, i] - mean(inno[, i]);
    inno[, i] <- inno[, i] / sd(inno[, i]);
    vola[, i] <- M@sigma.t;
    ics[i, ] <- M@fit$ics;
}
C <- cor(inno);

W <- matrix(NA, nrow=100*n, ncol=p);
sig2 <- matrix(NA, nrow=dim(W)[1], ncol=p);
# set the initial values
for (i in 1:p) {
    sig2[1, i] <- params[i, 1]/(1 - params[i, 3]);
    ## sig2[1, i] <- 0;
}
for (i in 1:dim(W)[1]) {
    ## eta <- rmvnorm(n=1, mean=rep(0, p), sigma=C);
    ## eta <- rmvt(n=1, sigma=C, df=params[, 4]);
    ## eta <- rt(n=dim(W)[2], df=params[, 4]);
    eta <- rnorm(n=dim(W)[2]);
    W[i, ] <- eta * sqrt(sig2[i,]);
    if (i < dim(W)[1]) {
        ## sig2[i+1, ] <- params[, 2] * W[i, ]^2 + params[, 1];
        sig2[i+1, ] <- params[, 2] * W[i, ]^2 +
            params[, 3] * sig2[i, ] + params[, 1];
    }
}
CX <- cov(X);
E <- eigen(CX);

CY <- cov(inno);
D <- eigen(CY);

CW <- cov(W);
F <- eigen(CW);

pdf("/tmp/GARCH-t_eigenvalues.pdf");
plot(1:p, E$values, type="p", pch=0,
     main="Spectra of FX and Simulated GARCH(1,1)-t",
     ylim=c(0, max(c(E$values, F$values))),
     xlab=expression(i),
     ylab=expression(lambda[(i)])
);
points(1:p, (F$values), pch=16, col="#FF0000");
legend("topright",
       legend=c(expression(cov(FX)), expression(cov(sim.^2))),
       ## col=c("#000000", "#0000FF", "#FF0000"),
       col=c("#000000", "#FF0000"),
       ## pch=c(0, 17, 16));
       pch=c(0, 16));
grid();
dev.off();

## pdf("/tmp/first_eigenvector.pdf");
## V <- sig.eig$vectors[, 1];
## V <- V * sign(V[which.max(abs(V))])
## plot(1:p, V, type="p", xaxt="n", xlab="i", ylab=expression(V[1](i)));
## axis(side=1, at=1:p, labels=names, las=2);
## dev.off();


pdf("/tmp/GARCH-t_eigenvectors.pdf", width=20, height=10);
par(mfrow=c(3,5));
## mse <- c(0, 0);
for (i in 1:p) {
    V <- E$vectors[, i];
    ## U <- D$vectors[, i];
    Q <- F$vectors[, i];
    
    ## if (sum(abs(V - U)) > sum(abs(V + U))) {
    ##     U <- -U;
    ## }
    if (sum(abs(V - Q)) > sum(abs(V + Q))) {
        Q <- -Q;
    }
    s <- sign(V[which.max(abs(V))]);

    ## mse[1] <- mse[1] + sum(abs(V * s - U * s));
    ## mse[2] <- mse[2] + sum(abs(V * s - Q * s));
    
    plot(1:p, V * s,
         main=sprintf("FX & sim. V[%d]", i),
         xlab="i", ylab=expression(V[i]),
         ylim=c(-1, 1), pch=0, cex=1.5,
         xaxt="n");
    axis(side=1, at=1:p, labels=names, las=2);
    
    points(1:p, Q * s, col="#FF0000", pch=16, cex=1.5);

    grid();
}
dev.off();
