# just a test

rm(list=ls());
library(R.matlab);
library(RMTstat)
# phi <- c(0, 0.5, 0.8160);
phi <- 0;
sig <- c(0.1, 0.2, 0.5);
# sig <- sqrt(1.0e-3);
# q <- 0.1;
q <- c(0.1, 0.2, 0.5, 1);
# TailIndices = matrix(nrow = length(q), ncol=length(sig));

## Values of T when sigma^2 = 1.0e-3;
# T <- c(2000, 2000, 1200, 600);
m = 1;
n = 1;
TWmean <- -1.2065335745820;
TWsd <- sqrt(1.607781034581);

X = seq(from=-4, to=4, length.out=2000);
Y = TWsd*dtw(TWsd * X + TWmean);
counter = 1;

for (n in 1:length(q)) {
    pdf(sprintf("./eigmax_q%.1f_dtw.pdf", q[n]));
    ##    par(mfrow=c(1, length(sig)));
    counter = 1;
    for (m in 1:length(sig)) {
        v <- sig[m]^2;
        data <- readMat(sprintf(paste("../matfys/data/sv/normal_ret/lognormal_vol/",
                                      "q%.1f/Eig-sig%.4f-phi%.4f.mat", sep=""),
                                q[n], sig[m], 0));
        eigmax = seq(length.out=dim(data$ev)[2]) * NA;
        for (k in 1:length(eigmax)) {
            eigmax[k] <- max(data$ev[, k])
        }
        eigmax = (eigmax - mean(eigmax)) / sd(eigmax);
        ## Par = WishartMaxPar(ndf=T[n], pdim=dim(data$ev)[1],
        ##     var=exp(2*v), beta=1);
        ## eigmax <- (eigmax - Par$centering) / Par$scaling;
        ## threshold = quantile(eigmax, 1 - 2.5e-2);
        ## E = sort(eigmax, decreasing=TRUE);
        ## ind = which(E > threshold);
        ## TailIndices[n, m] = 1/mean(log(E[ind] / threshold));

        epdf = density(eigmax);
        if (counter == 1) {
            plot(epdf$x, epdf$y, type="l", col=sprintf("#%06X", bitwShiftL(0xFF, 8*(m-1))),
                 ## main=sprintf("v=%.4f, tail index = %.4f", v, TailIndices[n, m]),
                 ## main=sprintf("v=%.4f, q=%.1f", v, q[n]),
                 main=sprintf("q=%.1f", q[n]),
                 xlab="X",
                 ylab="Densities");
                 ## col="blue", type="l");
        } else {
            lines(epdf$x, epdf$y, type="l", col=sprintf("#%06X", bitwShiftL(0xFF, 8*(m-1))));
        }
        counter = counter + 1;
        ## Y = dWishartMax(epdf$x, data$T, dim(data$ev)[1], exp(2*v));
    }
    lines(X, Y, col="#000000");
    dev.off();
}
## legend("topright", legend=c("v=0.01", "v=0.04", "v=0.25"));
## dev.off();





