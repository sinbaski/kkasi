rm(list=ls())
R <- read.table("path_stat.dat");
pdf("/tmp/right_eigenfunction.pdf");
plot(R[, 1], R[, 2], type="l",
     main=expression(r[xi]),
     xlab="x", ylab="y");
dev.off();

