rm(list=ls());
library(RMySQL);
library("fGarch");
library("MASS");
source("libxxie.r");

day2 = '2015-02-28';
day1 = '2010-01-01';

assetSet <- "SP500_components";

## database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
##     dbname='avanza', host=Sys.getenv("PB"));
database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host="localhost");
results = dbSendQuery(database, sprintf("select t1.symbol, t2.sector from
SP500_components as t1 join company_info as t2 on t1.symbol = t2.symbol where t2.sector='Materials';",
    assetSet));
H <- fetch(results, n=-1);
tables <- H[[1]];
sectors <- H[[2]];
dbClearResult(results);
n.stocks <- length(tables);
dbDisconnect(database);

X <- getInterpolatedReturns(day1, day2, "", tables, "_US");
n <- dim(X)[1];
p <- dim(X)[2];
res <- matrix(NA, nrow=n, ncol=p);
coef <- matrix(NA, nrow=p, ncol=3);
indices <- {};
for (i in 1:p) {
    ## M <- garch(x=X[, i], order=c(1, 1), trace=FALSE);
    ## coef[i, ] <- M$coef;
    ## res[, i] <- M$residuals;
    ## res[1, i] <- sign(rnorm(1));
    ## vol[, i] <- X[, i] / res[, i];

    M <- tryCatch({
         garchFit(~garch(1,1), data=X[, i], trace=FALSE);
    }, error = function(err) {
        FALSE;
    });
    if (typeof(M) != "logical") {
        coef[i, ] <- M@fit$params$params[c(2,3,5)];
        if (sum(coef[i, 2:3]) < 1) {
            res[, i] <- M@residuals;
            indices <- union(indices, i);
        }
    }

    ## M <- estimGARCH(0, 0.01, 0, X[, i]);
    ## coef[i, ] <- M$coef;
    ## res[, i] <- M$residus;
    ## vol[, i] <- X[11:n, i] / res[, i];
    ## print(c(names[i], coef[i, 2:3]));
}
res <- res[, indices];
coef <- coef[indices, ];
X <- X[, indices];
p <- length(indices);

C <- cor(res);
V <- apply(res, MARGIN=2, FUN=sd);
res <- res / matrix(rep(V, n), byrow=TRUE, nrow=n, ncol=p);
sigma <- X / res;

W <- matrix(NA, nrow=100*n, ncol=p);
sig2 <- matrix(NA, nrow=dim(W)[1], ncol=p);
# set the initial values
for (i in 1:p) {
    sig2[1, i] <- coef[i, 1]/(1 - coef[i, 3]);
    ## sig2[1, i] <- 0;
}
for (i in 1:dim(W)[1]) {
    eta <- mvrnorm(n=1, mu=rep(0, p), Sigma=C);
    W[i, ] <- eta * sqrt(sig2[i,]);
    if (i < dim(W)[1])
        sig2[i+1, ] <- coef[, 2] * W[i, ]^2 + coef[, 3] * sig2[i, ] + coef[, 1];
}

X2 <- X;
W2 <- W;

CX <- cov(X2);
E <- eigen(CX);
CW <- cov(W2);
F <- eigen(CW);

pdf("../papers/FX/Materials_eigenvalues.pdf");
plot(1:p, E$values, type="p", pch=17,
     xlab="i",
     ylab=expression(lambda[(i)]),
     main="Spectra of Real & Simulated Data", col="#000000"
);
points(1:p, (F$values), pch=16, col="#FF0000");

## ## points(1:p, (E1$values)/sum(E1$values), col="#FF0000", cex=2, pch=15);
## ## points(1:p, (F1$values)/sum(F1$values), col="#00FF00", cex=2, pch=16);
## ## points(1:p, (F1$values)/sum(F1$values), col="#00FF00", cex=2, pch=17);

legend("topright",
       legend=c(expression(cov("Materials")), expression(cov("Simulated"))),
       col=c("#000000", "#FF0000"),
       pch=c(17, 16), cex=2);
grid();
dev.off();

pdf("../papers/FX/Materials_eigenvectors1.pdf", width=10, height=10);
par(mfrow=c(4,3));
for (i in 1:12) {
    V <- E$vectors[, i];
    U <- F$vectors[, i];
    
    if (sum(abs(V - U)) > sum(abs(V + U))) {
        U <- -U;
    }
    s <- sign(V[which.max(abs(V))]);

    plot(1:p, V * s, main=sprintf("Materials & GARCH(1,1) V[%d]", i),
         xlab="i", ylab=expression(V[i]),
         ylim=c(-1, 1), pch=0,
         xaxt="n");
    axis(side=1, at=1:p, labels=tables[indices], las=2);
    
    points(1:p, U * s, main=sprintf("eigenvector[%d]", i),
           col="#FF0000", pch=16);
    grid();
}
dev.off();

pdf("../papers/FX/Materials_eigenvectors2.pdf", width=10, height=10);
par(mfrow=c(4,3));
for (i in 13:p) {
    V <- E$vectors[, i];
    U <- F$vectors[, i];
    
    if (sum(abs(V - U)) > sum(abs(V + U))) {
        U <- -U;
    }
    s <- sign(V[which.max(abs(V))]);

    plot(1:p, V * s, main=sprintf("Materials & GARCH(1,1) V[%d]", i),
         xlab="i", ylab=expression(V[i]),
         ylim=c(-1, 1), pch=0,
         xaxt="n");
    axis(side=1, at=1:p, labels=tables[indices], las=2);
    
    points(1:p, U * s, main=sprintf("eigenvector[%d]", i),
           col="#FF0000", pch=16);
    grid();
}
dev.off();
