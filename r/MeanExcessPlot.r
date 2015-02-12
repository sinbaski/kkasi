library(R.matlab);
library(RMTstat);

rm(list = ls());
phi <- c(0, 0.25, 0.5, 0.7070, 0.8160);
sig <- c(0.1, 0.2, 0.5);
# q <- c(0.1, 0.2, 0.5, 1);
q <- 0.1;

TailIndices = matrix(nrow = length(sig), ncol=length(phi));

pdf(sprintf("./eigmax_MeanExcess.pdf"));
# par(mfrow=c(length(sig), length(phi));
layout(matrix(1:15, 3, 5, byrow=TRUE));
for (m in 1:length(sig)) {
    for (n in 1 : length(phi)) {
        v = sig[m]^2 / (1 - phi[n]^2);
        data <- readMat(sprintf(paste("../matfys/data/sv/normal_ret/lognormal_vol/",
                                      "q%.1f/Eig-sig%.4f-phi%.4f.mat", sep=""),
                                q, sig[m], phi[n]));
        eigmax = seq(length.out=dim(data$ev)[2]) * NA;
        for (k in 1:length(eigmax)) {
            eigmax[k] <- max(data$ev[, k])
        }
        Par = WishartMaxPar(ndf=data$T, pdim=dim(data$ev)[1],
            var=exp(2*v), beta=1);
        eigmax = (eigmax - Par$centering) / Par$scaling;
        E = sort(eigmax, decreasing=FALSE);
        l = length(E);
        N = floor(l / 5);
        # thresholds = E[l - N : l - 1];
        meanExcess = rep(NA, N);

        for (i  in (l-N):(l-1)) {
            meanExcess[i-l+N+1] = mean(E[(i+1):l]  - E[i]);
        }
        plot(E[(l-N):(l-1)], meanExcess, type="p", pch=20,
             col="blue", xlab="", ylab="",
             main=sprintf("v=%.4f", v));
    }
}
dev.off();
