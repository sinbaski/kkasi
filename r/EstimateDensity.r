library(R.matlab);

X <- readMat("/tmp/DataSample.mat");
den <- density(X$sample);
## write.table(cbind(Den$x, Den$y), file="/tmp/Density.csv",
##             row.names=FALSE, col.names=FALSE);
density <- cbind(den$x, den$y);
writeMat("/tmp/DensityFunction.mat", density=density);

