rm(list=ls());
graphics.off();
library(RMySQL);
## library("fGarch");
require("rugarch");
require("mvtnorm");
source("libxxie.r");

day2 = '2015-02-28';
day1 = '2010-01-01';

assetSet <- "SP500_components";

## database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
##     dbname='avanza', host=Sys.getenv("PB"));
database = dbConnect(MySQL(), user='root', password='q1w2e3r4',
    dbname='avanza', host="localhost");
results = dbSendQuery(database, sprintf("select t1.symbol, t2.sector from
SP500_components as t1 join company_info as t2 on t1.symbol = t2.symbol where t2.sector='Utilities';",
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
coef <- matrix(NA, nrow=p, ncol=3);
inno <- matrix(NA, nrow=n, ncol=p);
indices <- {};
for (i in 1:p) {
    ## M <- garch(x=X[, i], order=c(1, 1), trace=FALSE);
    ## coef[i, ] <- M$coef;
    ## res[, i] <- M$residuals;
    ## res[1, i] <- sign(rnorm(1));
    ## vol[, i] <- X[, i] / res[, i];

    M <- tryCatch({
         garchFit(~garch(1,1),
                  data=X[, i],
                  trace=FALSE,
                  cond.dist="std",
                  shape=6,
                  include.shape=FALSE,
                  include.mean=TRUE,
                  include.delta=FALSE,
                  include.skew=FALSE
                  );
    }, error = function(err) {
        FALSE;
    });
    if (typeof(M) != "logical") {
        coef[i, ] <- M@fit$coef[-1];
        if (sum(coef[i, -1]) < 1) {
            indices <- union(indices, i);
            inno[, i] <- M@residuals / M@sigma.t;
            inno[, i] <- inno[, i] - mean(inno[, i]);
            inno[, i] <- inno[, i] / sd(inno[, i]);
        }
    }

    ## M <- estimGARCH(0, 0.01, 0, X[, i]);
    ## coef[i, ] <- M$coef;
    ## res[, i] <- M$residus;
    ## vol[, i] <- X[11:n, i] / res[, i];
    ## print(c(names[i], coef[i, 2:3]));
}
inno <- inno[, indices];
coef <- coef[indices, ];
X <- X[, indices];
p <- length(indices);

C <- cor(inno);

W <- matrix(NA, nrow=100*n, ncol=p);
sig2 <- matrix(NA, nrow=dim(W)[1], ncol=p);
# set the initial values
for (i in 1:p) {
    sig2[1, i] <- coef[i, 1]/(1 - coef[i, 3]);
    ## sig2[1, i] <- 0;
}
for (i in 1:dim(W)[1]) {
    eta <- rmvnorm(n=1, mean=rep(0, p), sigma=C);
    W[i, ] <- eta * sqrt(sig2[i,]);
    if (i < dim(W)[1])
        sig2[i+1, ] <- coef[, 2] * W[i, ]^2 + coef[, 3] * sig2[i, ] + coef[, 1];
}

CX <- cov(X);
E <- eigen(CX);
CW <- cov(W);
F <- eigen(CW);
CY <- cov(res);
D <- eigen(CY);

pdf("../papers/FX/Materials_eigenvalues.pdf");
plot(1:p, E$values/sum(E$values), type="p", pch=0,
     main="FX and GARCH(1,1) spectrum"
);
points(1:p, (D$values)/sum(D$values), col="#0000FF", pch=17);
points(1:p, (F$values)/sum(F$values), pch=16, col="#FF0000");
legend("topright",
       legend=c(expression(cov(real)), expression(cov(inno)), expression(cov(sim.))),
       col=c("#000000", "#0000FF", "#FF0000"),
       pch=c(0, 17, 16));
grid();
dev.off();

names <- tables[indices];
pdf("../papers/FX/Materials_eigenvectors1.pdf", width=10, height=10);
mse <- c(0, 0);
par(mfrow=c(5,5));
for (i in 1:25) {
    V <- E$vectors[, i];
    U <- D$vectors[, i];
    Q <- F$vectors[, i];
    
    if (sum(abs(V - U)) > sum(abs(V + U))) {
        U <- -U;
    }
    if (sum(abs(V - Q)) > sum(abs(V + Q))) {
        Q <- -Q;
    }
    s <- sign(V[which.max(abs(V))]);

    mse[1] <- mse[1] + sum(abs(V * s - U * s));
    mse[2] <- mse[2] + sum(abs(V * s - Q * s));
    
    plot(1:p, V * s, main=sprintf("FX & GARCH(1,1) V[%d]", i),
         xlab="i", ylab=expression(V[i]),
         ylim=c(-1, 1), pch=0,
         xaxt="n");
    axis(side=1, at=1:p, labels=names, las=2);
    
    ## points(1:p, U * s,
    ##        col="#0000FF", pch=17);
    points(1:p, Q * s, col="#FF0000", pch=16);

    grid();
}
dev.off();

pdf("../papers/FX/Materials_eigenvectors2.pdf", width=10, height=10);
par(mfrow=c(4,3));
for (i in 12:p) {
    V <- E$vectors[, i];
    U <- D$vectors[, i];
    Q <- F$vectors[, i];
    
    if (sum(abs(V - U)) > sum(abs(V + U))) {
        U <- -U;
    }
    if (sum(abs(V - Q)) > sum(abs(V + Q))) {
        Q <- -Q;
    }
    s <- sign(V[which.max(abs(V))]);

    mse[1] <- mse[1] + sum(abs(V * s - U * s));
    mse[2] <- mse[2] + sum(abs(V * s - Q * s));
    
    plot(1:p, V * s, main=sprintf("FX & GARCH(1,1) V[%d]", i),
         xlab="i", ylab=expression(V[i]),
         ylim=c(-1, 1), pch=0,
         xaxt="n");
    axis(side=1, at=1:p, labels=, las=2);
    
    points(1:p, U * s,
           col="#0000FF", pch=17);
##    points(1:p, Q * s, col="#FF0000", pch=16);

    grid();
}
dev.off();
