rm(list=ls());
library(R.matlab);
library(RMTstat)
                                        # phi <- c(0, 0.5, 0.8160);
phi <- 0;
## sig <- sqrt(c(0.01, 0.02, 0.03, 0.04, 0.25, 0.5, 1));
## sig <- c(0.1, 0.2, 0.5);
sig <- sqrt(seq(0.01, 0.02, by=0.002));
                                        # q <- 0.1;
q <- c(0.1, 0.2, 0.5, 1);
                                        # TailIndices = matrix(nrow = length(q), ncol=length(sig));

## Values of T when sigma^2 = 1.0e-3;
                                        # T <- c(2000, 2000, 1200, 600);
m = 1;
n = length(q);
TWmean <- -1.2065335745820;
TWsd <- sqrt(1.607781034581);

X = seq(from=0.5, to=4, length.out=2000);
Y = ptw(TWsd * X + TWmean, lower.tail=FALSE);
counter = 1;
colors = c("blue", 'cyan', "deepskyblue", "green", "darkgreen", 'red', 'darkorchid');
legends = rep("", length(sig));
## for (n in 1) {
pdf(sprintf("./eigmax_q%.1f_dtw_loglog_1.pdf", q[n]));
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

    F = ecdf(eigmax);
    legends[m] <- sprintf("v=%.3f", v);
    if (counter == 1) {
        plot(X, 1 - F(X), type="l",
             col=colors[m],
             pch=m,
             log="xy",
             main=sprintf("q=%.1f", q[n]),
             xlab="x",
             ylab="Tail Function");
        ## col="blue", type="l");
    } else {
        lines(X, 1 - F(X), type="l", pch=m, col=colors[m]);
    }
    counter = counter + 1;
    ## Y = dWishartMax(epdf$x, data$T, dim(data$ev)[1], exp(2*v));
}
lines(X, Y, col="black");
legend("topright", legend=c(legends, "TW"),
       lty=1, col=c(colors[1:length(sig)], "black"));
dev.off();
                                        # }

## dev.off();





