graphics.off();
H <- hillPlot(V[, 1]^2, prob=0.92, view=c(1.5, 2.5));
## H <- hillPlot(X[, 12]^2, prob=0.90);
labs <- seq(from=0, to=n/10, by=10);
ticks <- seq(from=1.5, to=3, by=0.2);
axis(side=1, at=labs);
axis(side=2, at=ticks);
abline(v=labs, h=ticks, lty=3);
H$alpha[80]
