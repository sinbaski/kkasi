rm(list=ls());
graphics.off();
## require("tseries");
require("fGarch");
require("mvtnorm");
require("stochvol")
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

# for t innovations
# params <- matrix(NA, nrow=p, ncol=4);
# for normal innovations
params <- matrix(NA, nrow=p, ncol=3);
inno <- matrix(NA, nrow=n, ncol=p);
H <- matrix(NA, nrow=n, ncol=p);
for (i in 1 : p) {
    ## result <- svsample(X[, i] - mean(X[, i]), priornu=c(3, 6));
    result <- svsample(X[, i] - mean(X[, i]));
    H[, i] <- apply(latent(result), MARGIN=2, FUN=mean);
    params[i, ] <- apply(para(result), MARGIN=2, FUN=mean);
    inno[, i] <- result$y / exp(H[, i]/2);
    inno[, i] <- inno[, i] - mean(inno[, i]);
}
C <- cor(inno);

Y <- matrix(NA, nrow=100*n, ncol=p);
H <- matrix(NA, nrow=100*n + 1, ncol=p);
Mu <- params[, 1];
Phi <- params[, 2];
Sigma <- params[, 3];
# normal innovations
H[1, ] <- rmvnorm(n=1, mean=Mu, sigma=diag(Sigma^2/(1 - Phi^2)));

for (t in 2 : dim(H)[1]) {
    H[t, ] <- rmvnorm(n=1, mean=Mu + Phi * (H[t-1, ] - Mu), sigma=diag(Sigma^2));
    Vol <- exp(H[t, ]/2);
    Y[t-1, ] <- rmvnorm(n=1, mean=rep(0, p), sigma=C);
    Y[t-1, ] <- Y[t-1, ] * Vol;
}

## params <- matrix(NA, nrow=p, ncol=2);
## inno <- matrix(NA, nrow=n, ncol=p);
## models <- vector("list", p);
## for (i in 1 : p) {
##     sv <- svsample(X[, i] - mean(X[, i]));
##     H <- apply(latent(sv), MARGIN=2, FUN=mean);
##     models[[i]] <- arima(H, order=c(1,2,1));
## ##     params[i, ] <- apply(para(sv), MARGIN=2, FUN=mean);
##     inno[, i] <- sv$y / exp(H/2);
##     inno[, i] <- inno[, i] - mean(inno[, i]);
## }
## C <- cor(inno);

## Y <- matrix(NA, nrow=100*n, ncol=p);
## for (i in 1 : p) {
##     Y[, i] <- arima.sim(models[[i]], dim(Y)[1]);
## }
## Y <- Y * rmvnorm(n=dim(Y)[1], mean=rep(0, p), sigma=C);

E <- eigen(cov(X));
D <- eigen(cov(inno));
F <- eigen(cov(Y));

pdf("/tmp/FX_sv_eigenvalues.pdf");
## plot(1:p, sig.eig$values, type="p", pch=17,
##      main="FX and GARCH(1,1) spectrum", col="#00FF00"
## );
plot(1:p, E$values/sum(E$values), type="p", pch=0,
     main="FX and SV spectrum",
     xlab=expression(i),
     ylab=expression(lambda[i]),
     cex=2,
     ylim=c(0, 1)
);
## points(1:p, (D$values)/sum(D$values), col="#0000FF", pch=17);
points(1:p, (F$values)/sum(F$values), pch=16, col="#FF0000", cex=2);

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

pdf("/tmp/FX_sv_eigenvectors.pdf", width=20, height=10);
par(mfrow=c(2,3));
mse <- c(0, 0);
for (i in 1:6) {
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
    
    plot(1:p, V * s, main=sprintf("FX & SV V[%d]", i),
         xlab="i", ylab=expression(V[i]),
         ylim=c(-1, 1), pch=0,
         cex=2,
         xaxt="n");
    axis(side=1, at=1:p, labels=names, las=2);
    
    ## points(1:p, U * s,
    ##        col="#0000FF", pch=17);
    points(1:p, Q * s, col="#FF0000", pch=16, cex=2);

    grid();
}
dev.off();
