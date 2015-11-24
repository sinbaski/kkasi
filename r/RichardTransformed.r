rm(list=ls());
source("libxxie.r");

day2 = '2015-02-28';
day1 = '2010-01-01';

assetSet <- "SP500_components";

R <- getInterpolatedReturns(day1, day2, assetSet);
R.trfm <- matrix(NA, nrow=dim(R)[1], ncol=dim(R)[2]);

for (i in 1 : length(stocks.included)) {
    ## Fn <- ecdf(R[, i]);
    ## U <- Fn(R[, i]);
    ## U <- U[U < 1];
    R.trfm[, i] <- -1/log(rank(R[,i])/(n.records+1));
}
# E <- eigen((n.records * p)^(-2) * t(R.trfm) %*% R.trfm);
E <- eigen((n.records * p)^(-2) * t(R) %*% R);

qtl <- function(u, i, a) {
    return(u^(2/(i*a)));
}

a <- 2.3;

## plot(lambda[2:p]/lambda[1:p-1], type="b", xlim=c(1, 20), ylim=c(0,1));

N = 50;
pdf(sprintf("../papers/Number1/EigenRatioSP500_log_%d_shown.pdf", N));

par(mar=c(4, 5, 2, 2));
plot((1:(N-1)), log(E$values[2:N]/E$values[1:(N-1)]), type="b", pch=19,
     xlim=c(1,N),
     ylim=(c(-2.0, 0)),
     xlab=expression(i),
     ylab=expression(log(lambda[(i+1)]/lambda[(i)])));
I <- 1:(N-1);
q = 0.99;
Q <- matrix(NA, ncol=3, nrow=N-1);
## Q[, 1] = q^(2/I);
Q[, 1] = qtl(0.99, I, a);
## Q[, 2] = (0.5)^(2/I);
Q[, 2] = qtl(0.5, I, a);
## Q[, 3] = (1-q)^(2/I);
Q[, 3] = qtl(0.01, I, a);

## expected <- I / (I+2);
lines((I), log(Q[, 2]), col="#0000FF", type="l", lwd=2);
lines((I), log(Q[, 3]), col="#FF0000", lwd=2);
lines((I), log(Q[, 1]), col="#00FF00", lwd=2);
dev.off();


