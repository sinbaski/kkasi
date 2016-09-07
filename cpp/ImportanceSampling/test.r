rm(list=ls())
a0 <- 1.0e-7;
a1 <- 0.11;
b1 <- 0.88;
a2 <- 1.0e-8;

norms <- rep(NA, 4000);

for (i in 1:4000) {
    Y <- rnorm(1)^2;
    A <- matrix(c(a1 * Y + b1, a2, Y, 0), byrow=TRUE, nrow=2, ncol=2);
    norms[i] <- norm(A, "I");
}

fun <- function(alpha) {
    return(mean(norms^alpha));
}

A <- seq(from=1, to=2, length.out=101);
X <- rep(NA, length(A));
for (i in 1:length(A)) {
    X[i] <- fun(A[i]);
}
pdf("A.pdf");
plot(A, X, type="l", xlab=expression(alpha), ylab=expression(E(group("||", A, "||"))^alpha));
dev.off();

