# just a test

rm(list=ls());
library(R.matlab);
library(RMTstat)
# phi <- c(0, 0.5, 0.8160);
# phi <- 0;
# sig <- c(0.1, 0.2, 0.5);
sig <- sqrt(1.0e-3);
# q <- 0.1;
q <- c(0.1, 0.2, 0.5, 1);
# TailIndices = matrix(nrow = length(q), ncol=length(sig));
T <- c(2000, 2000, 1200, 600);
m = 1;
n = 3;
pdf(sprintf("./eigmax_TW_dtw.pdf"));
par(mfrow=c(2, 2));
for (n in 1:length(q)) {
    for (m in 1:length(sig)) {
        v <- sig[m]^2;
        data <- readMat(sprintf(paste("../matfys/data/sv/normal_ret/lognormal_vol/",
                                      "q%.1f/Eig-sig%.4f-phi%.4f.mat", sep=""),
                                q[n], sig[m], 0));
        eigmax = seq(length.out=dim(data$ev)[2]) * NA;
        for (k in 1:length(eigmax)) {
            eigmax[k] <- max(data$ev[, k])
        }
        Par = WishartMaxPar(ndf=T[n], pdim=dim(data$ev)[1],
            var=exp(2*v), beta=1);
        eigmax <- (eigmax - Par$centering) / Par$scaling;
        ## threshold = quantile(eigmax, 1 - 2.5e-2);
        ## E = sort(eigmax, decreasing=TRUE);
        ## ind = which(E > threshold);
        ## TailIndices[n, m] = 1/mean(log(E[ind] / threshold));

        epdf = density(eigmax);
        plot(epdf$x, epdf$y,
             ## main=sprintf("v=%.4f, tail index = %.4f", v, TailIndices[n, m]),
             main=sprintf("v=%.4f, q=%.1f", v, q[n]),
             xlab="X",
             ylab="Densities",
             col="blue", type="l");
        ## Y = dWishartMax(epdf$x, data$T, dim(data$ev)[1], exp(2*v));
        Y = dtw(epdf$x);
        lines(epdf$x, Y, col="red");
    }
}
dev.off();





