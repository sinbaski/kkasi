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
N <- 8;

E <- cov(X);
Eig <- eigen(E);
C <-Eig$vectors[, 1:N];

Y <- X %*% C;

params <- matrix(0, nrow=N, ncol=3);
for (i in 1:N) {
    M <- garchFit(~garch(1,1),
                  data=Y[, i],
                  trace=FALSE,
                  ## cond.dist="std",
                  ## shape=4,
                  include.shape=FALSE,
                  include.mean=FALSE,
                  include.delta=FALSE,
                  include.skew=FALSE
                  );
    params[i, ] <- coef(M);
}

W <- matrix(NA, nrow=100*n, ncol=N);
sig2 <- matrix(NA, nrow=dim(W)[1], ncol=N);
# set the initial values
for (i in 1:N) {
    sig2[1, i] <- params[i, 1]/(1 - params[i, 3]);
    ## sig2[1, i] <- 0;
}
for (i in 1:dim(W)[1]) {
    eta <- rmvnorm(n=1, mean=rep(0, N), sigma=diag(rep(1, N)));
##    eta <- rmvt(n=1, sigma=C, df=4);
    W[i, ] <- eta * sqrt(sig2[i,]);
    if (i < dim(W)[1])
        sig2[i+1, ] <- params[, 2] * W[i, ]^2 + params[, 3] * sig2[i, ] + params[, 1];
}

D <- cov(W);
## Dig <- eigen(D);
    
## C <- matrix(NA, nrow=N, ncol=p);
## Res <- matrix(NA, nrow=n, ncol=p);
## for (i in 1:p) {
##     R <- X[, i];
##     model <- lm(R ~ Y[, 1:N]);
##     C[, i] <- coef(model)[1:N+1];
##     Res[, i] <- residuals(model);
## }
## Dep <- cov(Res);
F <- C %*% D %*% t(C);
Fig <- eigen(F);

pdf("/tmp/FX_OGARCH_eigenvalues.pdf", width=10, height=10);
plot(1:p, Eig$values, type="p", pch=0,
     xlab=expression(i),
     ylab=expression(lambda[(i)]),
     cex=2,
     main="FX and sim. spectrum", ylim=c(0, max(Eig$values[1], Fig$values[1]))
);
points(1:p, (Fig$values), pch=16, col="#FF0000", cex=2);

## ## points(1:p, (E1$values)/sum(E1$values), col="#FF0000", cex=2, pch=15);
## ## points(1:p, (F1$values)/sum(F1$values), col="#00FF00", cex=2, pch=16);
## ## points(1:p, (F1$values)/sum(F1$values), col="#00FF00", cex=2, pch=17);

legend("topright",
##       legend=c(expression(sigma[i] * sigma[j]), expression(cov(W)), expression(cov(X))),
       legend=c(expression(cov(FX)), expression(cov(sim.))),
       col=c("#000000", "#FF0000"),
       pch=c(0, 16));
grid();
dev.off();

pdf("/tmp/FX_OGARCH_eigenvectors.pdf", width=20, height=10);
par(mfrow=c(2,4));
## mse <- c(0, 0);
for (i in 1:8) {
    V <- Eig$vectors[, i];
    Q <- Fig$vectors[, i];
    
    if (sum(abs(V - Q)) > sum(abs(V + Q))) {
        Q <- -Q;
    }
    s <- sign(V[which.max(abs(V))]);

    ## mse[1] <- mse[1] + sum(abs(V * s - U * s));
    ## mse[2] <- mse[2] + sum(abs(V * s - Q * s));
    
    plot(1:p, V * s, main=sprintf("FX & sim. V[%d]", i),
         xlab="i", ylab=expression(V[i]),
         ylim=c(-1, 1), pch=0, cex=2,
         xaxt="n");
    axis(side=1, at=1:p, labels=names, las=2);
    
    ## points(1:p, U * s,
    ##        col="#0000FF", pch=17);
    points(1:p, Q * s, col="#FF0000", pch=16, cex=2);

    grid();
}
dev.off();
