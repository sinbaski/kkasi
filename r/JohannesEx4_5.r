rm(list=ls());

qMyDist <- function(U, a) {
    X <- rep(NA, length(U));
    I1 <- which(U <= 1/4)
    U1 <- U[I1];
    I2 <- which(U > 1/4 & U <= 3/4);
    U2 <- U[I2];
    I3 <- which(U > 3/4);
    U3 <- U[I3];
    X[I1] <- -(4^(a+1) * U1)^(-1/a);
    X[I2] <- U2/a - (1 + 1/a)/4;
    X[I3] <- 4^(1 + 1/a) * (1 - U3)^(1/a);
    X[I3] <- 1/X[I3];
    return(X);
}

rMyDist <- function(n, a) {
    U <- runif(n);
    return (qMyDist(U, a));
}


p <- 200;
n <- 1000;

N <- 1000;
a <- 0.6;
lambda <- matrix(rep(NA, N*2), nrow=N, ncol=2);
for (i in (1:1000)) {
    ## Z <- matrix(rMyDist((p+1)*(n+1), a), nrow=p+1, ncol=n+1);
    ## X <- Z[2:(p+1), 2:(n+1)] + Z[2:(p+1), 1:n] - 2 * (Z[1:p, 2:(n+1)] - Z[1:p, 1:n]);
    X <- matrix(rMyDist(p*n, a), nrow=p, ncol=n);
    C <- X %*% t(X);
    E <- eigen(C, only.values=TRUE);
    lambda[i, ] <- E$values[1:2];
}
W <- (lambda[, 1] - lambda[, 2])/lambda[, 1];
Fn <- ecdf(W);
X <- seq(0, 1, length.out=400);
Y <- Fn(X);
pdf("lambda_1_2_Gaps_iid.pdf")
plot((X), (Y), type="b", xlab="x",
     ylab=expression(P((lambda[(1)] - lambda[(2)])/lambda[(1)] < x)),
     main="Distribution Function of the self-normalized gap");
## I <- X <= 3/4;
## Y1 <- rep(NA, length(X));
## Y1[I] <- 1 - (1 - X[I])^(a/2);
## Y1[!I] <- 1;
## lines(X[I], Y1[I], col="#FF0000", type="l", lwd=2);
## lines(X[!I], Y1[!I], col="#FF0000", type="l", lwd=2);
lines(X, 1 - (1 - X)^(a/2), col="#FF0000", type="l", lwd=2);
grid(nx=20);

dev.off();


