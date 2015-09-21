rm(list=ls());

source("liblxb353.r");

p <- 200;
n <- 1000;

## The eigenvector of the largest eigenvalue of a Wishart matrix is delocalozed.
                                        # X <- matrix(rMyDist(p * n, a=0.8), nrow=p, ncol=n);
X <- matrix(rnorm(p * n), nrow=p, ncol=n);
C <- X %*% t(X) / n;
E <- eigen(C);
T <- E$vectors[, 1];

pdf("LargestEigenvector-NormalData.pdf");
plot(1:length(T), T, type="b",
     main="N(0,1) data",
     xlab="Indices of components",
     ylab="Size of components");
dev.off();





