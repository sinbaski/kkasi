## Simulates eigenvalues of fixed-dimensional sample covariance matrix
## of correlated infinite MA processes
rm(list=ls());

source("liblxb353.r");

p <- 200;
n <- 1000;

mu <- (sqrt(n-1) + sqrt(p))^2;
sigma <- (sqrt(n-1) + sqrt(p)) * (1/sqrt(n-1) + 1/sqrt(p))^(1/3);

## Simulate using myDist distribution
N <- 2000;
a <- 1.6;
a.np <- qMyDist(1 - 1/(n*p), a) * 1.6;
lambda <- rep(NA, N);
for (i in (1001:2000)) {
    X <- matrix(rMyDist(p*n, a), nrow=p, ncol=n);
    C <- X %*% t(X);
    E <- eigen(C, only.values=TRUE);
    lambda[i] <- E$values[1];
}
W <- lambda[1:N]/a.np^2;
W1 <- (lambda - mu)/sigma;

Fn <- ecdf(W);
X <- exp(seq(min(log(W)), max(log(W)), length.out=200));
Y <- Fn(X);

pdf("MyDist-Frechet.pdf")
plot(log(X), (Y), type="p", xlab=expression(log(x)),
     ylab=expression(P(lambda[(1)] < x)),
     main="Sample Distribution Function and Frechet");
Y.f <- pfrechet(X, shape=a/2);
lines(log(X), Y.f, col="#FF0000");
grid(nx=20);
explanations <- c(
    "sample",
    "Frechet"
    );
legend("topleft", legend=explanations,
       lty=c(1, 1), lwd=c(1, 1),
       col=c("#000000", "#FF0000"));
dev.off();



## Simulate using the 3-point distribution
## N <- 2000;
## lambda <- rep(NA, N);
## for (i in (2:1000)) {
##     X <- matrix(r3Point(p*n, 1/6, sqrt(3)), nrow=p, ncol=n);
##     C <- X %*% t(X);
##     E <- eigen(C, only.values=TRUE);
##     lambda[i] <- E$values[1];
## }


## Simulate using N(0, 1)
## N <- 2000;
## lambda <- rep(NA, N);
## for (i in (1001:2000)) {
##     X <- matrix(rnorm(p*n), nrow=p, ncol=n);
##     C <- X %*% t(X);
##     E <- eigen(C, only.values=TRUE);
##     lambda[i] <- E$values[1];
## }

# W <- (lambda[1:1000] - mu)/sigma;

## den <- density(W);

## pdf(sprintf("3point-TW.pdf"));
## plot(den$x, den$y, type="l", xlab="x", ylab="density",
##      main="Sample Density function and Tracy-Widom");
## Y <- dtw(den$x);
## lines(den$x, Y, col="#FF0000");
## grid(nx=20);
## explanations <- c(
##     "sample",
##     "Tracy-Widom"
##     );
## legend("topleft", legend=explanations,
##        lty=c(1, 1), lwd=c(1, 1),
##        col=c("#000000", "#FF0000"));
## dev.off();

## Simulate with myDist: uniform in the center and symmetric pareto on the tails
## alpha <- 1.6;
## for (i in (1:100)) {
##     Anp <- (n*p)^(1/alpha);
##     data <- rMyDist(p*n, alpha);
##     X <- matrix(data, nrow=p, ncol=n);
##     C <- X %*% t(X);
##     E <- eigen(C/Anp^2, only.values=TRUE);
##     D <- sort(rowSums(X^2/Anp^2), decreasing=TRUE);
##     Zupp <- sort((data/Anp)^2, decreasing=TRUE)[1:p];
##     save(E, D, Zupp, file=sprintf("Eigen_D_Nbr%d_alpha%.1f_n%d_p%d.dat",
##                    i, alpha, n, p), compress="gzip");
## }




## for (n in seq(1000, 2000, by=200)) {
##     p <- 5*floor(log(n));
##     a.n <- qt(1-1/n, alpha);
##     fold <- 200;

##     ## iid student dist. sample _fold_ lambda.1 values
##     lambda.1 <- rep(NA, fold);
##     for (i in 1:fold) {
##         R <- matrix(rt(p*n, alpha), nrow=p, ncol=n)  / a.n;
##         C <-  R %*% t(R);
##         E <- eigen(C);
##         A <- log(E$values);
##         A <- (A - mean(A))/sd(A);
##         lambda.1[i] <- exp(A[1]);
##     }
##     save("lambda.1", file=sprintf("lambda.1.alpha%.1f_n%d.dat", alpha, n),
##          compress="gzip");
## }
