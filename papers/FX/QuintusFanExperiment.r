rm(list=ls());
require("fGarch");
source("libxxie.r");


h <- 0.49;
## set.seed(19);
Gn <- QuintosFanRollingDist(h, 8000, 1000);

n <- 2000;
p <- 1000;
## set.seed(9);
## Z <- matrix(rmydist(n=n*p, 0.8), nrow=n, ncol=p);
## Z <- matrix(NA, nrow=n, ncol=p);
## spec <- garchSpec(model=list(alpha=0.11, beta=0.88));
## for (i in 1:p) {
##     Z[, i] <- garchSim(spec, n=n);
## }

Z <- matrix(NA, nrow=n, ncol=p);
index <- 2.5;
for (i in 1:p) {
    n1 <- floor(n/2);
    Z[1:n1, i] <- rt(n=n1, index);
    Z[(n1 + 1):n, i] <- rt(n=n - n1, index - 0.5);
}


Q <- apply(Z, MARGIN=2,
           FUN=function(S) {
               max(QuintosFanRollingTest(
                   S, k=0.03,
                   h=h, variant="iid"))
           });
S <- Gn(Q);
sum(S > 0.99);
