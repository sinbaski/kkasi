rm(list=ls())
a0 <- 1.0e-7;
a1 <- 0.6;
b1 <- 0.005;
a2 <- 0.001;

norms <- rep(NA, 4000);

for (i in 1:4000) {
    Y <- rnorm(1)^2;
    A <- matrix(c(a1 * Y + b1, a2, Y, 0), byrow=TRUE, nrow=2, ncol=2);
    norms[i] <- norm(A, "I");
}

fun <- function(alpha) {
    return(mean(norms^alpha));
}

X <- seq(from=0, to=1.5, length.out=101);
Y <- rep(NA, length(X));
for (i in 1:length(X)) {
    Y[i] <- log(fun(X[i]));
}
## pdf("A.pdf");
plot(X, Y, type="l", xlab=expression(alpha), ylab=expression(E(group("||", A, "||")^alpha)));
grid();

## dev.off();

