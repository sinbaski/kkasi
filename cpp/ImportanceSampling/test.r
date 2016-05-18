rm(list=ls())

## for (alpha in c(1.5, 2.0, 2.5, 3.0, 3.5)) {
##     X <- read.table(sprintf("Lambda_%2.1f.txt", alpha))$V1;
##     cat(sprintf("% 2.2f    % 6.4f    %6.4f\n", alpha, mean(X), sd(X)));
## }

data <- read.table("dat/garch21-2.dat", skip=3, comment.char=c('#'));
lo <- loess(data$V2 ~ data$V1);
pdf("/tmp/GARCH21-2.pdf")
plot(data$V1, data$V2, type="l", xlab=expression(alpha),
     ylab=expression(Lambda(alpha)),
     main=expression(
         sigma[t+1]^2 == 10^{-7} + 0.11*R[t]^2
         +0.09*R[t-1]^2 + 0.79*sigma[t]^2
     ));
lines(data$V1, predict(lo), col="#00FF00");
grid();
dev.off();
## X <- read.table("Lambda_2.0.txt")$V1;
## Y <- read.table("Lambda_1.5.txt")$V1;
## h <- shapiro.test(X);
## h$p.value
## dens <- density(X);
## plot(dens$x, dens$y);
## qqnorm(X);

