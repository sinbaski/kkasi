rm(list=ls())
R <- read.table("eta.dat");
## pdf("../../papers/Jeffrey1/SP500_r.pdf");
plot(R[, 1], R[, 2], type="p",
     main=expression(pi),
     xlab="x", ylab="y");
## dev.off();

