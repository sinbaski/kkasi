rm(list=ls())

## for (alpha in c(1.5, 2.0, 2.5, 3.0, 3.5)) {
##     X <- read.table(sprintf("Lambda_%2.1f.txt", alpha))$V1;
##     cat(sprintf("% 2.2f    % 6.4f    %6.4f\n", alpha, mean(X), sd(X)));
## }

data <- read.table("garch21.dat", skip=3, comment.char=c('#'));
pdf("/tmp/GARCH21.pdf")
plot(data$V1, data$V2, type="l", xlab=expression(alpha),
     ylab=expression(Lambda(alpha)),
     main=expression(
         sigma[t+1]^2 == 10^{-7} + 0.11*R[t]^2
         +10^{-8}*R[t-1]^2 + 0.88*sigma[t]^2
     ));
grid();
dev.off();
## X <- read.table("Lambda_2.0.txt")$V1;
## Y <- read.table("Lambda_1.5.txt")$V1;
## h <- shapiro.test(X);
## h$p.value
## dens <- density(X);
## plot(dens$x, dens$y);
## qqnorm(X);

