rm(list=ls());
graphics.off();
source("libxxie.r");

omega <- 1.0e-7;
a <- c(0.11, 1.0e-8);
b <- 0.88;

n <- 500000;
V <- matrix(NA, nrow=2, ncol=n);
V[, 1] <- c(omega, 0);

A <- matrix(NA, 2, 2);
B <- c(omega, 0);

for (i in 2:n) {
    Z2 <- rnorm(1)^2;
    A[1, 1] <- a[1] * Z2 + b;
    A[1, 2] <- a[2];
    A[2, 1] <- Z2;
    A[2, 2] <- 0;
    V[, i] <- A %*% V[, i-1] + B;
}

Ka <- hillPlot(V[2, 10000:n]);
plot(Ka$K, Ka$alpha, ylim=c(1.5, 2.5), xlim=c(0, 4000), type="l");
abline(h=seq(from=1000, to=4000, by=200), v=0:3, col="#FF0000");

## h <- hillEstimate(X[2, ]);
