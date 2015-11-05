rm(list=ls());

alpha <- 1.6;
p <- 200;
n <- 1000;
k <- sqrt(p);

fold <- 1000;
lambda <- matrix(NA, nrow=p, ncol=fold);
Ds <- matrix(NA, nrow=p, ncol=fold);
Zupps <- matrix(NA, nrow=p, ncol=fold);
for (i in (1:fold)) {
    load(file=sprintf("dat/Eigen_D_Nbr%d_alpha%.1f_n%d_p%d.dat",
             i, alpha, n, p));
    lambda[, i] <- E$values;
    Ds[, i] <- D;
    Zupps[, i] <- Zupp;
}

## S <- matrix(rep(NA, p*2), nrow=2, ncol=p);
## for (i in 1:p) {
##     ## S[i] <- max(abs(Ds[i,]/lambda[i,] - 1));
##     S[1, i] <- max(abs(Ds[i,] - lambda[i,]));
##     S[2, i] <- median(abs(Ds[i,] - lambda[i,]));
## }
## pdf("lambda_comparison1.pdf");
## plot(1:p, S[1,], xlim=c(1,50), xlab="i", ylab="", type="l", lty=1, col="#FF0000");
## lines(1:p, S[2,], type="l", lty=4, col="#0000FF");
## explanations <- c(expression(a[n*p]^{-2} * max[j] * bgroup("|", D[(i)]^{(j)} - lambda[(i)]^{(j)}, "|")),
##                   expression(a[n*p]^{-2} * median[j] * bgroup("|", D[(i)]^{(j)} - lambda[(i)]^{(j)}, "|")));
## legend("topright", legend=explanations, lty=c(1, 4), cex=1.5);

## dev.off();

pdf("../papers/Number1/lambda_comparison.pdf");
errors <- lambda[1,] - Ds[1,];
a <- quantile(errors, prob=0.02);
b <- quantile(errors, prob=0.98);
den <- density(errors[which(errors > a & errors < b)]);
denfun <- approxfun(den$x, den$y*0.96);
plot(den$x, den$y, type="l", xlab="", ylab="Distribution density",
     col="#FF0000", lwd=2, lty=1);
x0 <- median(errors);
y0 <- denfun(x0);
y <- seq(0, y0, by=y0/100);
lines(rep(x0, 101), y, col="#FF0000", lty=2);

errors <-  lambda[1,] - Zupps[1,];
a <- quantile(errors, prob=0.02);
b <- quantile(errors, prob=0.98);
den <- density(errors[which(errors > a & errors < b)]);
denfun <- approxfun(den$x, den$y*0.96);
x0 <- median(errors);
y0 <- denfun(x0);
y <- seq(0, y0, by=y0/100);
lines(den$x, den$y, col="#0000FF", lwd=2, lty=4)
lines(rep(x0, 101), y, col="#0000FF", lty=2);

grid(nx=20);

explanations <- c(
    expression(a[np]^{-2}*(lambda[(1)] - D[(1)]^{"" %->% ""})),
    "median",
    expression(a[np]^{-2}*(lambda[(1)] - Z[(1)]^2)),
    "median"
    );
legend("topright", legend=explanations,
       col=c("#FF0000", "#FF0000", "#0000FF", "#0000FF"),
       lty=c(1,2,4,2), lwd=c(2,1,2,1), cex=2);

dev.off();

## pdf("lambda_minus_Z_squared.pdf");
## subject <- expression(a[np]^{-2}*(Z[(i)]^2 - lambda[(i)]));
## ## hist(lambda - Zupps, breaks=200, xlab=subject, main="Histogram of approximation error",
## ##      xlim=c(-1.0e-3, 2.0e-3), ylim=c(0, 18000));
## errors <- Zupps - lambda;
## a <- quantile(errors, prob=0.01);
## b <- quantile(errors, prob=0.99);
## den <- density(errors[which(errors > a & errors < b)]);
## plot(den$x, den$y, type="l", xlab=subject, ylab="Distribution density");
## grid(nx=20);
## dev.off();
