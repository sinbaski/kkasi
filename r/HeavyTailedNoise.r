rm(list=ls());
## library(fields);
source("libxxie.r");

#Simulate samples described by a heavy-tailed SV model.
n <- 100000;
p <- 20;


## Z <- runif(1.0e+5);
## Y <- (Z^(-1/3));
## X <- matrix(rt(n=n*p, df=5), nrow=n, ncol=p);
X <- runif(n*p)^(-1/3);
X <- X * sign(runif(n*p, min=-0.5, 0.5));
X <- matrix(X, n, p);
## estimateTailIndices(X);

Y <- matrix(NA, n, p);
eta <- matrix(rnorm(n=(p+1)*n), nrow=n, ncol=p+1);
for (i in 2:(p+1)) {
    Y[, i-1] <- X[, i-1] * exp(eta[, i-1] + eta[,i]);
}

tail.indices <- matrix(NA, p, p);
for (i in (1:p)) {
    for (j in (i:p)) {
        tail.indices[i,j] <- pikandsEstimate(Y[,i] * Y[,j], prob=0.97);
        tail.indices[j,i] <- tail.indices[i,j]
    }
}

M = max(tail.indices);
m <- min(tail.indices);

colors <- gray((tail.indices - m)/(M - m));
## my.colors <- colorRampPalette(c("#00008B", "#8B0000"));
## colors <- my.colors(p*(p+1)/2);

pdf("/tmp/ParetoNoise_AbsnormVolatility_Hill.pdf")
plot(1, 1, type="n", xlim=c(1,p+2), ylim=c(1,p));
for (i in 1:p) {
    for (j in i:p) {
        points(x=i, y=j, pch=19, cex=3, col=colors[(j-1)*p+i]);
        ## points(x=i, y=j, pch=19, cex=3,
        ##        col=colors[floor((tail.indices[i,j]-m)/(M-m)*length(colors))+1]);
    }
}
## axis(side=1, at=1:p, labels=names);
## axis(side=2, at=1:p, labels=names);
## image.plot(legend.only=TRUE, zlim=c(m, M), col=sort(colors))
image.plot(x=seq(1,p+2), y=1:p, legend.only=TRUE, zlim=c(m, M), col=sort(colors));
dev.off();

