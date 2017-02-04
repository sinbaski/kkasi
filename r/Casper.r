rm(list=ls());
graphics.off();
## require("tseries");
## require("fGarch");
## require("mvtnorm");
source("libxxie.r");

## currencies <- c(
##     "AUD_SEK_Rates",
##     "CAD_SEK_Rates",
##     "CHF_SEK_Rates",

##     ## "CNY_SEK_Rates",
##     "CZK_SEK_Rates",
##     "DKK_SEK_Rates",

##     "EUR_SEK_Rates",
##     "GBP_SEK_Rates",
##     "HKD_SEK_Rates",

##     "HUF_SEK_Rates",
##     "JPY_SEK_Rates",
##     "KRW_SEK_Rates",

##     ## "MAD_SEK_Rates",
##     ## "MXN_SEK_Rates",
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
##     "CHF",

##     ## "CNY",
##     "CZK",
##     "DKK",

##     "EUR",
##     "GBP",
##     "HKD",

##     "HUF",
##     "JPY",
##     "KRW",

##     ## "MAD",
##     ## "MXN",
##     "NOK",

##     "NZD",
##     ##"SAR",
##     "SGD",
##     "USD"
##     );

## X <- getAssetReturns("2010-01-03", "2016-04-01", currencies, 1,
##                      "rate", "localhost");

tables <- c(
    "ADS_DE",
    "CBK_DE",
    "DBK_DE",
    "DTE_DE",
    "EOAN_DE",
    "SIE_DE",
    "VOW3_DE"
);

X <- getInterpolatedReturns("2008-01-01", "2014-06-30", "", tables, "");

n <- dim(X)[1];
p <- dim(X)[2];
## N <- 9;
## X <- X - matrix(rep(apply(X, FUN=mean, MARGIN=2), n), nrow=n, ncol=p, byrow=TRUE);

E <- eigen(cov(X));
C <-E$vectors[, 1:N];
Y <- X %*% C;

inno <- matrix(NA, nrow=n, ncol=N);
vola <- matrix(NA, nrow=n, ncol=N);
params <- matrix(0, nrow=N, ncol=3);
params.cvar <- matrix(NA, nrow=N, ncol=4);
## ARCH
## params <- matrix(0, nrow=p, ncol=2);
ics <- matrix(NA, nrow=p, ncol=4);
for (i in 1:N) {
    M <- garchFit(~garch(1,1),
                  data=Y[, i],
                  trace=FALSE,
                  ## cond.dist="std",
                  cond.dist="norm",
                  ## shape=4,
                  include.shape=FALSE,
                  include.mean=FALSE,
                  include.delta=FALSE,
                  include.skew=FALSE
                  );
    params[i, ] <- coef(M);
    params.cvar[i, ] <- as.vector(M@fit$cvar[-1, -1]);
    inno[, i] <- M@residuals / M@sigma.t;
    ## inno[, i] <- inno[, i] - mean(inno[, i]);
    ## inno[, i] <- inno[, i] / sd(inno[, i]);
    vola[, i] <- M@sigma.t;
    ics[i, ] <- M@fit$ics;
}

hill.params <- c(2.911810e-02,	9.699083e-01,
                 4.269052e-02,	9.552375e-01,
                 5.276025e-02,	9.436635e-01,
                 4.226736e-02,	9.561051e-01,
                 8.174961e-02,	9.081800e-01,
                 2.996618e-02,	9.690467e-01,
                 1.737636e-02,	9.823525e-01,
                 3.770294e-02,	9.606773e-01,
                 3.096452e-02,	9.680478e-01
                 );

hill.params <- matrix(hill.params, nrow=length(hill.params)/2,
                      ncol=2, byrow=TRUE);
params[1:dim(hill.params)[1], c(2,3)] = hill.params;

## N <- dim(hill.params)[1];
## params <- params[1:N, ];

V <- matrix(NA, nrow=n*100, ncol=N);
sig2 <- matrix(1, nrow=dim(V)[1], ncol=N);
for (i in 1:N) {
    sig2[1, i] <- params[i, 1]/(1 - params[i, 3]);
}
for (i in 1:dim(V)[1]) {
    ## eta <- rt(n=N, df=params[, 4]);
    ## eta <- rmvt(n=1, sigma=C, df=df);
    eta <- rmvnorm(n=1, mean=rep(0, N),
                   sigma=diag(rep(1, N)));
    V[i, ] <- eta * sqrt(sig2[i,]);
    if (i < dim(V)[1])
        ## sig2[i+1, ] <- fixed.garch[, 2] * V[i, ]^2 + fixed.garch[, 1];
        sig2[i+1, ] <- params[, 2] * V[i, ]^2 +
            params[, 3] * sig2[i, ] + params[, 1];
}

W <- V %*% t(C);
F <- eigen(cov(W));
# F <- eigen(C %*% cov(V) %*% t(C));
pdf("/tmp/FX_eigenvalues.pdf");
plot(1:p, E$values/sum(E$values), type="p", pch=0,
     main="Spectra of FX and Simulated GARCH(1,1) Series",
     ylim=c(0, 1),
     xlab=expression(i),
     ylab=expression(lambda[(i)])
);
points(1:p, (F$values)/sum(F$values), pch=16, col="#FF0000");
legend("topright",
       legend=c(expression(cov(FX)), "Simulated GARCH"),
       ## col=c("#000000", "#0000FF", "#FF0000"),
       col=c("#000000", "#FF0000"),
       ## pch=c(0, 17, 16));
       pch=c(0, 16));
grid();
dev.off();

pdf("/tmp/FX_eigenvectors.pdf", width=20, height=10);
par(mfrow=c(2,5));
## mse <- c(0, 0);
for (i in 1:N) {
    E.evec <- E$vectors[, i];
    ## U <- D$vectors[, i];
    F.evec <- F$vectors[, i];
    
    ## if (sum(abs(E.evec - U)) > sum(abs(E.evec + U))) {
    ##     U <- -U;
    ## }
    if (sum(abs(E.evec - F.evec)) > sum(abs(E.evec + F.evec))) {
        F.evec <- -F.evec;
    }
    s <- sign(E.evec[which.max(abs(E.evec))]);

    ## mse[1] <- mse[1] + sum(abs(E.evec * s - U * s));
    ## mse[2] <- mse[2] + sum(abs(E.evec * s - F.evec * s));
    
    plot(1:p, E.evec * s,
         main=sprintf("FX & Simulated E.evec[%d]", i),
         xlab="i", ylab=expression(E.evec[i]),
         ylim=c(-1, 1), pch=0, cex=1.5,
         xaxt="n");
    axis(side=1, at=1:p, labels=names, las=2);
    
    points(1:p, F.evec * s, col="#FF0000", pch=16, cex=1.5);

    grid();
}
dev.off();
