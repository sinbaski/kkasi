rm(list=ls());

alpha <- 1.6;
p <- 200;
n <- 1000;
k <- sqrt(p);

fold <- 1000;
lambda <- matrix(NA, nrow=p, ncol=fold);
Ds <- matrix(NA, nrow=p, ncol=fold);
Zupps <- matrix(NA, nrow=p, ncol=fold);
Ys <- matrix(NA, nrow=p, ncol=fold);
for (i in (1:fold)) {
    load(file=sprintf("dat/iid/Eigen_D_Nbr%d_alpha%.1f_n%d_p%d.dat",
             i, alpha, n, p));
    lambda[, i] <- E$values;
    Ds[, i] <- D;
    Zupps[, i] <- Zupp;
    Ys[,i] <- sort(Y, decreasing=TRUE);
}

## pdf("/tmp/comparison2.pdf");
## A <- Zupps - Ys;
## errors <- rep(NA, nrow=1, ncol=fold);
## for (i in (1:fold)) {
##     errors[i] <- max(abs(A[,i]));
## }
## a <- quantile(errors, prob=0.02);
## b <- quantile(errors, prob=0.98);
## den <- density(errors[which(errors > a & errors < b)]);
## denfun <- approxfun(den$x, den$y*0.96);
## plot(den$x, den$y, col="#FF0000", type="l", lwd=2, lty=1,
##      xlab=expression(a[np]^{-2}*max(bgroup("[", Z[(i)]^2 - Y[(i)], "]"), i<=p)),
##      ylab="Distribution density");
## dev.off();



pdf("/tmp/comparison2.pdf");
errors1 <-  lambda[1,] - Ds[1,];
## errors1 <- rep(NA, nrow=1, ncol=fold);
## for (i in (1:fold)) {
##     errors1[i] <- max(abs(A[,i]));
## }
a <- quantile(errors1, prob=0.02);
b <- quantile(errors1, prob=0.98);
den <- density(errors1[which(errors1 > a & errors1 < b)]);
## den <- density(errors1);
denfun <- approxfun(den$x, den$y*0.96);
x0 <- median(errors1);
y0 <- denfun(x0);
y <- seq(0, y0, by=y0/100);
plot(den$x, den$y, col="#FF0000", type="l", lwd=2, lty=1,
     xlab="", ylab="Smoothed Histogram");
# lines(rep(x0, 101), y, col="#FF0000", lty=2);

errors2 <-  lambda[1,] - Zupps[1,];
## errors2 <- rep(NA, nrow=1, ncol=fold);
## for (i in (1:fold)) {
##     errors2[i] <- max(abs(A[,i]));
## }
a <- quantile(errors2, prob=0.02);
b <- quantile(errors2, prob=0.98);
den <- density(errors2[which(errors2 > a & errors2 < b)]);
denfun <- approxfun(den$x, den$y*0.96);
## den <- density(errors2);
## denfun <- approxfun(den$x, den$y);
x0 <- median(errors2);
y0 <- denfun(x0);
y <- seq(0, y0, by=y0/100);
lines(den$x, den$y, col="#0000FF", type="l", lwd=4, lty=4)
## plot(den$x, den$y, col="#0000FF", type="l", lwd=2)
## lines(rep(x0, 101), y, col="#0000FF", lty=2);

## A <-  lambda[,] - Ys[,];
## errors3 <- rep(NA, nrow=1, ncol=fold);
## for (i in (1:fold)) {
##     errors3[i] <- max(abs(A[,i]));
## }
## a <- quantile(errors3, prob=0.02);
## b <- quantile(errors3, prob=0.98);
## den <- density(errors3);
## denfun <- approxfun(den$x, den$y);
## x0 <- median(errors3);
## y0 <- denfun(x0);
## y <- seq(0, y0, by=y0/100);
## lines(den$x, den$y, col="#00FF00", type="l", lwd=2)


explanations <- c(
    expression(a[np]^{-2}*group("(", lambda[(1)] - D[(1)]^{"" %->% ""}, ")")),
    expression(a[np]^{-2}*group("(", lambda[(1)] - Z[(1)]^2, ")"))
    );
## explanations <- c(
##     expression(a[np]^{-2}*max(bgroup("|", lambda[(i)] - D[(i)]^{"" %->% ""}, "|"), i<=p)),
##     expression(a[np]^{-2}*max(bgroup("|", lambda[(i)] - Z[(i)]^2, "|"), i<=p))
##     );
## legend("topright", legend=explanations,
##        col=c("#FF0000", "#FF0000", "#0000FF", "#0000FF"),
##        lty=c(1,2,4,2), lwd=c(2,1,2,1), cex=2);
legend("topright", legend=explanations,
       col=c("#FF0000", "#0000FF"),
       lty=c(1,4), lwd=c(2,2), cex=2);

dev.off();


## pdf("../papers/Number1/lambda_comparison.pdf");
## errors <- lambda[1,] - Ds[1,];
## a <- quantile(errors, prob=0.02);
## b <- quantile(errors, prob=0.98);
## den <- density(errors[which(errors > a & errors < b)]);
## denfun <- approxfun(den$x, den$y*0.96);
## plot(den$x, den$y, type="l", xlab="", ylab="Distribution density",
##      col="#FF0000", lwd=2, lty=1);
## x0 <- median(errors);
## y0 <- denfun(x0);
## y <- seq(0, y0, by=y0/100);
## lines(rep(x0, 101), y, col="#FF0000", lty=2);

## errors <-  lambda[1,] - Zupps[1,];
## a <- quantile(errors, prob=0.02);
## b <- quantile(errors, prob=0.98);
## den <- density(errors[which(errors > a & errors < b)]);
## denfun <- approxfun(den$x, den$y*0.96);
## x0 <- median(errors);
## y0 <- denfun(x0);
## y <- seq(0, y0, by=y0/100);
## lines(den$x, den$y, col="#0000FF", lwd=2, lty=4)
## lines(rep(x0, 101), y, col="#0000FF", lty=2);

## grid(nx=20);

## explanations <- c(
##     expression(a[np]^{-2}*(lambda[(1)] - D[(1)]^{"" %->% ""})),
##     "median",
##     expression(a[np]^{-2}*(lambda[(1)] - Z[(1)]^2)),
##     "median"
##     );
## legend("topright", legend=explanations,
##        col=c("#FF0000", "#FF0000", "#0000FF", "#0000FF"),
##        lty=c(1,2,4,2), lwd=c(2,1,2,1), cex=2);

## dev.off();

