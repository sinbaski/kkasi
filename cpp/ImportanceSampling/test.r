rm(list=ls())

## for (alpha in c(1.5, 2.0, 2.5, 3.0, 3.5)) {
##     X <- read.table(sprintf("Lambda_%2.1f.txt", alpha))$V1;
##     cat(sprintf("% 2.2f    % 6.4f    %6.4f\n", alpha, mean(X), sd(X)));
## }

X <- read.table("inf_norms_n300_pow2.txt")$V1;
## X <- read.table("Lambda_2.0.txt")$V1;
## Y <- read.table("Lambda_1.5.txt")$V1;
## h <- shapiro.test(X);
## h$p.value
## dens <- density(X);
## plot(dens$x, dens$y);
## qqnorm(X);

