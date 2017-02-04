rm(list=ls());
require("mvtnorm");
source("libxxie.r");

set.seed(0);
C <- matrix(NA, nrow=15, ncol=2);
C[, 1] <- (1:15)*0.1;
C[, 2] <- rep(1, 15);
## C <- rep(1, 15);
N <- dim(C)[1] + 1;
n <- 1500;
indices <- matrix(NA, nrow=dim(C)[1], ncol=2);
## indices <- rep(NA, length(C));

params <- matrix(rep(c(1.366814e-06, 0.03310482, 0.9358017), N),
                 nrow=N, ncol=3, byrow=TRUE);

V <- matrix(NA, nrow=n, ncol=N);
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

## X <- matrix(NA, nrow=n, ncol=N-1);
for (i in 1:dim(C)[1]) {
    ## X[, i] <- C[i] * V[, i] + V[, N];
    X <- C[i, 1] * V[, i] + V[, N];
    indices[i, 1] <- hillEstimate(X^2);
    X <- C[i, 2] * V[, i] + V[, N];
    indices[i, 2] <- hillEstimate(X^2);
}
## for (i in 1:dim(X)[2]) {
##     indices[i] <- hillEstimate(X[, i]^2);
## }

graphics.off();
pdf("/tmp/simulated_indices.pdf");
plot(1:dim(indices)[1], indices[, 1], ylab="Indices", pch=15);
points(1:dim(indices)[1], indices[, 2], pch=1, col="#FF0000");
abline(v=1:dim(indices)[1], h=seq(from=2.5, by=0.5, to=5), lty=3);
legend("topleft", legend=c("scale from 0.1 to 1.5", "scale=1"),
       col=c("#000000", "#FF0000"), pch=c(15, 1));
dev.off();
