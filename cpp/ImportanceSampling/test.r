rm(list=ls())

X <- read.table("data.txt")$V1;
h <- shapiro.test(X);
h$p.value
dens <- density(X);
plot(dens$x, dens$y);
qqnorm(X);

