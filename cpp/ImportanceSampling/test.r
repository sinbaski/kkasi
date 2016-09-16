rm(list=ls())
a0 <- 1.0e-7;
a1 <- 0.11;
b1 <- 0.88;
a2 <- 1.0e-8;

norms <- rep(NA, 4000);

## title <- sprintf("a1 = %.3f, a2 = %.3f, b1 = %.3f", a1, a2, b1);
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

pdf("/tmp/upper_bound.pdf");
plot(X, Y, type="l",
     main=expression(sigma[t]^2 == 0.11 * R[t-1]^2 + 10^{-8} * R[t-2]^2 + 0.88 * sigma[t-1]^2 + 10^{-7}),
     xlab=expression(theta),
     ylab=expression(E(group("||", A, "||")^theta)));
grid();
dev.off();

