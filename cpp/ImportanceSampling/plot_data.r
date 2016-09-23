rm(list = ls());
graphics.off();
## data <- read.table(file="Lambda1.txt");
## pdf("/tmp/GARCH21_1.pdf")
## plot(data$V1, data$V2, type="l", xlab=expression(theta),
##      ylab=expression(Lambda(theta)),
##      main=expression(sigma[t]^2 == 0.6 * R[t-1]^2 + 10^{-3} * R[t-2]^2 + 0.005 * sigma[t-1]^2 + 10^{-7}));
## lines(data$V1, data$V2 + data$V3, col="#FF0000");
## lines(data$V1, data$V2 - data$V3, col="#0000FF");
## grid();
## dev.off();

old <- read.table(file="./lambda0.txt");
new <- read.table(file="./lambda0_new.txt");

pdf("old_and_new_algo.pdf")
plot(old$V1, old$V2,
     type="l", xlab=expression(theta),
     ylab=expression(Lambda(theta)),
     main=expression(sigma[t]^2 == 0.11 * R[t-1]^2 + 0.88 * sigma[t-1]^2 + 10^{-7})
     );
lines(new$V1, new$V2, col="#00FF00");
legend("topleft",
       legend=c("old", "new"),
       col=c("#000000", "#00FF00"),
       lwd=c(1, 1),
       );
grid();
dev.off();