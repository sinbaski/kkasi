rm(list=ls());
source("libxxie.r");
## source("liblxb353.r");

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


pdf("Eigenvalues_Cor_vs_Cov_Matrix.pdf")
plot(log10(E$values), log10(E.adjusted$values), ylim=c(0,5), xlim=c(-2.5, 2.5),
     xlab=expression(log[10](lambda^{cov})), ylab=expression(log[10](lambda^{corr})));
grid(nx=10);
dev.off();


## pdf("Eigenvalues_Cor_n_Cov_Matrix.pdf")
## plot(1:p, (E.adjusted$values), col="#FF0000", type="b",
##      ylim=c(0, 1.0E4),
##      xlab="i", ylab=expression(lambda[(i)]));
## points(1:p, (E$values), type="b", col="#000000");
## grid(nx=10);
## legend("topright", legend=c("Corr. Matrix", "Cov. Matrix"),
##        pch=c(1, 1), col=c("#FF0000", "#000000"));
## dev.off();


for (i in 1:6) {
    pdf(sprintf("V%d_CovarianceMatrix.pdf", i));
    plot(1:p, E$vectors[,i], xlab="k",
         ylab="kth comp",
         main=sprintf("Eigenvector of the %dth largest eigenvalue", i));
    grid(nx=10);
    dev.off();

    pdf(sprintf("V%d_CorrelationMatrix.pdf", i));
    plot(1:p, E.adjusted$vectors[,2], xlab="k",
         ylab="kth comp",
         main=sprintf("Eigenvector of the %dth largest eigenvalue", i));
    grid(nx=10);
    dev.off();
}

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

