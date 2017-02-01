rm(list=ls());
graphics.off();
source("libxxie.r");

tables <- c(
    "ADS_DE",
    "CBK_DE",
    "DBK_DE",
    "DTE_DE",
    "EOAN_DE",
    "LHA_DE",
    "SIE_DE",
    "VOW3_DE"
);

X <- getInterpolatedReturns("2008-01-01", "2014-06-30", "", tables, "");

n <- dim(X)[1];
p <- dim(X)[2];

prob <- 0.02;
t0 <- 0.4;
Fn <- QuintosFanRollingDist(t0, n.paths=2000, n.steps=1000);

H <- matrix(NA, nrow=p, ncol=p);
for (i in 1:(p-1)) {
    for (j in (i+1):p) {
        h <- QuintosFanRollingTest(c(X[, i]^2, X[, j]^2), prob, t0);
        ## the p-value of the test result
        H[i, j] <- 1 - Fn(max(h));
        H[j, i] <- H[i, j];
    }
}
diag(H) <- 1;
