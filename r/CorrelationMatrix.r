rm(list=ls());
source("libxxie.r");

day2 = '2015-02-28';
day1 = '2010-01-01';

assetSet <- "SP500_components";

R <- getInterpolatedReturns(day1, day2, assetSet);

n <- dim(R)[1];
p <- dim(R)[2];

### The covariance matrix
E <- eigen(t(R) %*% R);

### The correlation matrix
R.adjusted <- R;
for (i in 1:p) {
    sigma <- sd(R.adjusted[,i]);
    R.adjusted[,i] <- (R.adjusted[,i] - mean(R.adjusted[,i]))/sigma;
}
E.adjusted <- eigen(t(R.adjusted) %*% R.adjusted);
E.adjusted$values <- abs(E.adjusted$values);


pdf("Eigenvalues_Cor_n_Cov_Matrix.pdf")
plot(1:p, log10(E.adjusted$values), col="#FF0000", type="b",
     xlab="i", ylab=expression(log[10](lambda[(i)])));
points(1:p, log10(E$values), type="b", col="#000000");
grid(nx=10);
legend("bottomleft", legend=c("Correlation Matrix", "Covariance Matrix"),
       pch=c(1, 1), col=c("#FF0000", "#000000"));
dev.off();

## pdf("V1_CovarianceMatrix.pdf");
## plot(1:p, E$vectors[,1], xlab="k", ylim=c(-0.1, -0.01),
##      ylab=expression(V[k1]),
##      main="Eigenvector of the largest eigenvalue");
## grid(nx=10);
## dev.off();

## pdf("V1_CorrelationMatrix.pdf");
## plot(1:p, E.adjusted$vectors[,1], xlab="k", ylim=c(-0.1, -0.01),
##      ylab=expression(V[k1]),
##      main="Eigenvector of the largest eigenvalue");
## grid(nx=10);
## dev.off();

## pdf("Eigenvalues_of_CovarianceMatrix_SP500.pdf");
## plot(1:(p-1), log10(E$values[2:p]/E$values[1:(p-1)]), type="b",
##      xlim=c(1,50), ylim=c(-1.25, 0),
##      xlab="i", ylab=expression(log[10](lambda[(i+1)]/lambda[(i)])),
##      main="Eigenvalue Ratios of SP500 Covariance Matrix");
## dev.off();

## pdf("LargestComponent_of_Eigenvectors_of_CovarianceMatrix_SP500.pdf");
## plot(1:p, apply(abs(E$vectors), 1, "max"),
##      main="Component of largest abs. value of cov. matrix's eigenvectors",
##      xlab="i", ylab="comp. of largest abs. value");
## dev.off();


## pdf("Eigenvalues_of_CorrelationMatrix_SP500.pdf");
## plot(1:(p-1), log10(E.adjusted$values[2:p]/E.adjusted$values[1:(p-1)]), type="b",
##      xlim=c(1,50), ylim=c(-1.25, 0),
##      xlab="i", ylab=expression(log[10](lambda[(i+1)]/lambda[(i)])),
##      main="Eigenvalue Ratios of SP500 Correlation Matrix");
## dev.off();

## pdf("LargestComponent_of_Eigenvectors_of_CorrelationMatrix_SP500.pdf");
## plot(1:p, apply(abs(E.adjusted$vectors), 1, "max"),
##      main="Component of largest abs. value of corr. matrix's eigenvectors",
##      xlab="i", ylab="comp. of largest abs. value");
## dev.off();

