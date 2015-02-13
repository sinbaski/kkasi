library(R.matlab);

data <- readMat(sprintf(paste("../matfys/data/sv/normal_ret/lognormal_vol/",
                              "q%.1f/Eig-sig%.4f-phi%.4f.mat", sep=""),
                        0.1, 0.1, 0));
eigmax = seq(length.out=dim(data$ev)[2]) * NA;
for (k in 1:length(eigmax)) {
    eigmax[k] <- max(data$ev[, k])
}
Par = WishartMaxPar(ndf=data$T, pdim=dim(data$ev)[1],
    var=exp(2*v), beta=1);
eigmax = (eigmax - Par$centering) / Par$scaling;
E = sort(eigmax, decreasing=FALSE);
N = floor(length(E) / 5);
thresholds = tail(E, N);
meanExcess = rep(NA, N);

for (i  in 1:N-1) {
    meanExcess[i] = mean(tail(E, N-i )  - thresholds[i]);
}
plot(thresholds, meanExcess, type="p");

