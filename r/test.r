# just a test

library(R.matlab);
library(RMTstat)
phi <- c(0, 0.5, 0.8160);
# phi <- 0;
sig <- c(0.1, 0.2, 0.5);
q <- 0.1;
# sig <- 0.5;
pdf(sprintf("eigmax_TW_q%.2f.pdf", q));
par(mfrow=c(3, 3));
for (m in 1:length(sig)) {
    for (n in 1:length(phi)) {

        v <- sig[m]^2/(1 - phi[n]^2);
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

        epdf = density(eigmax);
        plot(epdf$x, epdf$y, main=sprintf("v=%.4f", v),
             xlab="X",
             ylab="Densities",
             col="blue", type="l");
        # X = seq(from=min(eigmax), to=max(eigmax), length.out=400);
        ## Y = dWishartMax(epdf$x, ndf=data$T, pdim=dim(data$ev)[1],
        ##     var=exp(v), beta=1);
        Y = dtw(epdf$x);
        lines(epdf$x, Y, col="red");
        
        ## eigmax1 <- rWishartMax(length(eigmax), ndf=data$T, pdim=dim(data$ev)[1],
        ##                        var=exp(v/2), beta=1);
        ## qqplot(eigmax, eigmax1, xlab=sprintf("v=%.4f", v), ylab="", col="blue", type="l");
    }
}
dev.off();





