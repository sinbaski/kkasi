rm(list=ls());
library(fields);
source("~/work/r/libxxie.r");

#Simulate samples described by a heavy-tailed SV model.
n <- 4000;
p <- 20;

X <- matrix(rt(n=n*p, df=3), nrow=p, ncol=n);
eta <- matrix(rnorm(n=(p+1)*n), nrow=p+1, ncol=n);
for (i in 2:(p+1)) {
    X[i-1,] <- X[i-1,] * exp(eta[i-1,] + eta[i]);
}

tail.indices <- array(NA, dim=c(p,p,2));
for (i in (1:p)) {
    for (j in (i:p)) {
        tail.indices[i,j,] <- hillEstimate(X[i,] * X[j,], probs=0.97);
        tail.indices[j,i,] <- tail.indices[i,j,]
    }
}

## lower.tail <- tail.indices[,,1];
## upper.tail <- tail.indices[,,2];

A <- tail.indices[,,2]

M <- max(A[which(!is.na(A), arr.ind=TRUE)]);
m <- min(A[which(!is.na(A), arr.ind=TRUE)]);

colors <- gray((A - m)/(M - m));
## my.colors <- colorRampPalette(c("#00008B", "#8B0000"));
## colors <- my.colors(p*(p+1)/2);

pdf("/tmp/HeavyTailedNoise_Hill.pdf")
plot(1, 1, type="n", xlim=c(1,p+2), ylim=c(1,p));
for (i in 1:p) {
    for (j in i:p) {
        points(x=i, y=j, pch=19, cex=3, col=colors[(i-1)*p+j]);
        ## points(x=i, y=j, pch=19, cex=3,
        ##        col=colors[floor((A[i,j]-m)/(M-m)*length(colors))+1]);
    }
}
## axis(side=1, at=1:p, labels=names);
## axis(side=2, at=1:p, labels=names);
## image.plot(legend.only=TRUE, zlim=c(m, M), col=sort(colors))
image.plot(x=seq(1,p+2), y=1:p, legend.only=TRUE, zlim=c(m, M), col=sort(colors));
dev.off();

