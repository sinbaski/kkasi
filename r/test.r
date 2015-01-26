# just a test

library(R.matlab);
library(RMTstat)
Matlab$startServer();
matlab <- Matlab();
open(matlab);
phi <- c(0, 0.5, 0.8160);
sig <- c(0.1, 0.2, 0.5);
q <- 0.1;

for (m in 1:length(sig)) {
    for (n in 1:length(phi)) {
        v <- sig[m]^2/(1 - phi[n]^2);
        data <- readMat(sprintf(paste("../matfys/data/sv/normal_ret/lognormal_vol/",
            "q%.1f/Eig-sig%.4f-phi%.4f.mat", sep=""), q, sig[m], phi[n]));
        dmsn = dim(data$R);
        eigmax = seq(length.out=dim(data$ev[2])) * NA;
        for (k in 1:length(eigmax)) {
            eigmax[k] = max(data$ev[1:dmsn[1]])
        }
        setVariable(data$ev);
        ev <- matrix(data$ev, nrow=1);
        X <- seq(min(ev), max(ev), length.out=400);
        Ye <- ecdf(ev)(X)
        Y <- pWishartMax(X, ndf=data$T, pdim=dim(data$ev)[1], var=exp(v/2),
            beta=1, lower.tail=TRUE, log.p=FALSE);
        pdf(sprintf("eigmax_dist_v%.3f.pdf", v));
        plot(X, Ye, col="blue");
        lines(X, Y, col="green");
        dev.off();
    }
}
close(matlab);
rm(list = ls())




