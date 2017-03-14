rm(list=ls())
R <- read.table("SP500_r.dat");
pdf("../../papers/Jeffrey1/SP500_r.pdf");
plot(R[, 1], R[, 2], type="l",
     main="right eigenfunction",
     xlab="x", ylab=expression(r[xi](x)));
dev.off();

