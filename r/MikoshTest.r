rm(list=ls());

alpha <- 1.6;
p <- 200;
n <- 1000;
k <- sqrt(p);

fold <- 100;
lambda <- matrix(NA, nrow=p, ncol=fold);
Ds <- matrix(NA, nrow=p, ncol=fold);
Zupps <- matrix(NA, nrow=p, ncol=fold);
for (i in (1:fold)) {
    load(file=sprintf("Eigen_D_Nbr%d_alpha%.1f_n%d_p%d.dat",
             i, alpha, n, p));
    lambda[, i] <- E$values;
    Ds[, i] <- D;
    Zupps[, i] <- Zupp;
}

pdf("lambda_comparison.pdf");
# subject <- expression(a[np]^{-2}*(D[(i)] - lambda[(i)]));
## hist(lambda - Zupps, breaks=200, xlab=subject, main="Histogram of approximation error",
##      xlim=c(-1.0e-3, 2.0e-3), ylim=c(0, 18000));
errors <- Zupps - lambda;
a <- quantile(errors, prob=0.01);
b <- quantile(errors, prob=0.99);
den <- density(errors[which(errors > a & errors < b)]);
denfun <- approxfun(den$x, den$y);
plot(den$x, den$y, type="l", xlab="", ylab="Distribution density",
     col="#0000FF", ylim=c(0, 26000), lwd=2, lty=1);
x0 <- median(errors);
y0 <- denfun(x0);
y <- seq(0, y0, by=y0/100);
lines(rep(x0, 101), y, col="#0000FF", lty=2);

errors <- Ds - lambda;
a <- quantile(errors, prob=0.01);
b <- quantile(errors, prob=0.99);
den <- density(errors[which(errors > a & errors < b)]);
denfun <- approxfun(den$x, den$y);
x0 <- median(errors);
y0 <- denfun(x0);
y <- seq(0, y0, by=y0/100);
lines(den$x, den$y, col="#FF0000", lwd=2, lty=4)
lines(rep(x0, 101), y, col="#FF0000", lty=2);

grid(nx=20);

explanations <- c(
    expression(a[np]^{-2}*(Z[(i)]^2 - lambda[(i)])),
    "median",
    expression(a[np]^{-2}*(D[(i)] - lambda[(i)])),
    "median"
    );
legend("topleft", legend=explanations,
       col=c("#0000FF", "#0000FF", "#FF0000", "#FF0000"),
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
