rm(list=ls());
library(FMStable)

alpha <- 0.8;
ns <- seq(2200, 5400, by=800);
supremum <- rep(NA, length(ns));

c <- 1;
for (n in ns) {
    p <- 5*floor(log(n));
    load(sprintf("dat/P-5xlogN_Dist-t/lambda.1.alpha%.1f_n%d.dat", alpha, n));
    X <- sort(lambda.1);
    supremum[c] <- max(c(abs(F[-p] - seq(1, p-1) / p), abs(F[-1] - seq(1, p-1) / p)));
    c <- c + 1;
}



## Let the H matrix be exponentially decaying along each row.
## A <- 1:p;
## C <- A;
## M <- 2*p;
## for (i in 1 : p)

