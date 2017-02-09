rm(list=ls())
R <- read.table("right_eigenfunction.dat");
plot(R[, 1], R[, 2], type="l",
     main="right eigenfunction",
     xlab="x", ylab=expression(r[alpha](x)));

