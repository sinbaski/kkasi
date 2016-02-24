rm(list=ls());
library(abind);

p <- 10;
n <- 10000;
## number of instances
M <- 200;

# eigenvalues divided by the tr(XX')
test1 <- matrix(NA, M, p);

# ratio of successive eigenvalues
test2 <- matrix(NA, M, p-1);

# differene of the largest and the smallest divided by tr(XX')
## test3 <- rep(NA, M);

test3 <- matrix(NA, M, 11);

psi <- c(1, 1);

for (iterator in 18) {
    n <- 2^iterator;
    sprintf("n=%d", n);
    for (i in 1:M) {
        Y <- matrix(rexp((p+1)*n, 1), p+1, n);
        Sigma <- matrix(NA, p, n);
        for (j in 1:p) {
            Sigma[j,] <- exp(psi %*% Y[j:(j+1),]);
        }
        X <- Sigma * matrix(rnorm(p*n), p, n);
        C <- (X %*% t(X)/n^2);
        tr <- sum(diag(C));
        E <- eigen(C, only.values=TRUE);
        ## test1[i, ] <- E$values/tr;
        ## test2[i, ] <- E$values[1:(p-1)] / E$values[2:p];
        test3[i, iterator-13] <- (E$values[1] - E$values[p])/tr;
    }
}
