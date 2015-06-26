## Simulates eigenvalues of fixed-dimensional sample covariance matrix
## of correlated infinite MA processes
rm(list=ls());
alpha <- 0.8;
for (n in seq(2200, 5400, by=800)) {
    p <- 5*floor(log(n));
    a.n <- qt(1-1/n, alpha);
    fold <- 200;

    ## iid student dist. sample _fold_ lambda.1 values
    lambda.1 <- rep(NA, fold);
    for (i in 1:fold) {
        R <- matrix(rt(p*n, alpha), nrow=p, ncol=n)  / a.n;
        C <-  R %*% t(R);
        E <- eigen(C);
        lambda.1[i] <- E$values[1];
    }
    save("lambda.1", file=sprintf("lambda.1.alpha%.1f_n%d.dat", alpha, n),
         compress="gzip");
}
