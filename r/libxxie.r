library(RMySQL);
library(abind);
## library(rmgarch);

###
# A matrix of return values derived from price data on common days
# of all stocks in the given tables.
### 
getAssetReturns <- function(day1, day2, tables, lag,
                            col.name, host) {
    R <- getAssetPrices(day1, day2, tables, lag, col.name, host);
    return(diff(log(R)));
}

getAssetPrices <- function(day1, day2, tables, lag,
                           col.name, host)
{
    database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
        dbname='avanza', host=host);
    days <- vector('character');
    for (i in 1:length(tables)) {
        results <- dbSendQuery(
            database,
            sprintf("select distinct day from %s where
day between '%s' and '%s' order by day;", tables[i], day1, day2));
        D = fetch(results, n=-1)[[1]];
        dbClearResult(results);
        if (length(days) > 0) {
            days <- intersect(days, D);
        } else {
            days <- D;
        }
    }
    n.days = length(days);
    R = matrix(nrow=ceiling(n.days/lag), ncol=length(tables));
    str = sprintf("'%s'", days[1]);
    for (i in 2 : n.days) {
        str = sprintf("%s, '%s'", str, days[i]);
    }

    for (i in 1:length(tables)) {
        results <- dbSendQuery(database, sprintf("select %s from
%s where day in (%s) order by day;", col.name, tables[i], str));
        prices <- fetch(results, n=-1)[[1]];
        dbClearResult(results);
        I <- rev(seq(from=length(prices), to=1, by=-lag));
        R[,i] <- prices[I]);
    }
    dbDisconnect(database);
    return (R);
}

### Given a vector of observations, infer the innovations
inferInnovations <- function(X) {
    p <- length(X);
    acov <- acf(X, type="covariance", plot=FALSE)$acf;
    acov <- append(acov, rep(0, length(X) - length(acov) + 1));
    InnAlgo <- innovations.algorithm(acov);
    inno <- rep(NA, p);
    inno[1] <- X[1];
    for (i in 1:p-1) {
        inno[i+1] <- X[i+1] - sum(InnAlgo$thetas[[p]][1:i]*rev(inno[1:i]));
    }
    return(inno);
}

getInterpolatedReturns <- function(day1, day2, assetSet) {
    database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
        dbname='avanza', host="localhost");
    results = dbSendQuery(database, sprintf("select symbol from %s;",
        assetSet));
    tables <- fetch(results, n=-1)[[1]];
    dbClearResult(results);
    n.stocks <- length(tables);
    suffix <- "";
    if (assetSet == "SP500_components") {
        suffix <- "_US";
    } else if (assetSet == "DAX_components") {
        suffix <- "_DE";
    }


    results <- dbSendQuery(
        database,
        sprintf("select count(*) from Calendar where weekday(day) <= 4
and day between '%s' and '%s';", day1, day2));
    n.records <- fetch(results, n=1)[1, 1];
    dbClearResult(results);

    to.include <- rep(TRUE, n.stocks);
    for (i in 1 : n.stocks) {
        tables[i] <- gsub("[.]", "_", tables[i]);
        tables[i] <- gsub("-", "_series_", tables[i]);
        tables[i] <- paste(tables[i], suffix, sep="");

        results <- dbSendQuery(
            database,
            sprintf("select count(*) from %s;",
                    tables[i])
            );
        num.records <- fetch(results)[1,1];
        if (num.records == 0) {
            to.include[i] <- FALSE;
            next;
        }
        
        results <- dbSendQuery(
            database,
            sprintf("select min(day) from %s;",
                    tables[i])
            );
        d1 <- fetch(results)[1, 1];
        dbClearResult(results);

        results <- dbSendQuery(
            database,
            sprintf("select max(day) from %s;", tables[i])
            );
        d2 <- fetch(results)[1,1];
        dbClearResult(results);

        if (d1 > day1 || d2 < day2) {
            to.include[i] <- FALSE;
            next;
        }

    }
    p <- sum(to.include);
    prices <- matrix(NA, nrow=n.records, ncol=p);

    results <- dbSendQuery(
        database,
        sprintf("select day from Calendar where weekday(day) <= 4
and day between '%s' and '%s';", day1, day2)
        );
    days <- fetch(results, n=-1)[[1]];
    dbClearResult(results);

    stocks.included <- which(to.include);
    for (i in 1:length(stocks.included)) {
        results <- dbSendQuery(
            database,
            sprintf("select day, closing from %s
where day between '%s' and '%s'",
                    tables[stocks.included[i]],
                    day1, day2)
            );
        A <- fetch(results, n=-1);
        dbClearResult(results);
        I <- days %in% A[[1]];
        prices[which(I), i] <- A[[2]];

        if (is.na(prices[1, i])) {
            ## get the last trading day before the period
            results <- dbSendQuery(
                database,
                sprintf("select max(day) from %s where day <= '%s';",
                        tables[stocks.included[i]], day1)
                );
            d <- fetch(results)[1, 1];
            dbClearResult(results);

            ## get the closing price of the last trading day
            results <- dbSendQuery(
                database,
                sprintf("select closing from %s where day = '%s';",
                        tables[stocks.included[i]], d)
                )
            last.price <- fetch(results)[1, 1];
            dbClearResult(results);

            ## fill in the price up to the first available
            j <- 1;
            while (is.na(prices[j, i])) {
                prices[j, i] <- last.price;
                I[j] <- TRUE;
                j <- j + 1;
            }
        }
        for (j in which(!I)) {
            prices[j, i] <- prices[j-1, i];
        }
    }
    dbDisconnect(database);
    R <- diff(log(prices));
}

estimateTailIndices <- function(ret, prob=0.97) {
    if (class(ret) == "matrix") {
        tailIndices <- rep(NA, dim(ret)[2]);
        for (i in 1:dim(ret)[2]) {
            X = ret[, i];
            b <- quantile(X, probs=prob);
            tailIndices[i] <- 1/mean(log(X[which(X > b)]/b));
        }
    } else if (class(ret) == "numeric") {
        b <- quantile(ret, probs=prob);
        tailIndices <- 1/mean(log(ret[which(ret > b)]/b));
    }
    return(tailIndices);
}

cor.conf.ntvl <- function (prob, N) {
    z <- qnorm(prob, mean=0, sd=1/sqrt(N-3));
    x <- exp(2*z);
    return ((x - 1)/(x + 1));
}

computeCovCorr <- function(data, max.lag=60) {
    lagged.cov <- array(dim=c(p,p,1));
    lagged.cor <- array(dim=c(p,p,1));
    r <- cor.conf.ntvl(c(0.01, 0.99), n);
    h <- 1;
    sigma <- apply(data, 2, "sd");
    sigma.inv <- sigma %o% sigma;
    while (h <= max.lag) {
        A <- matrix(ncol=p, nrow=p);
        for (i in 1:p) {
            for (j in 1:p) {
                A[i, j] <- cov(data[1:(n-h), i], data[(1+h):n, j]);
            }
        }
        B <- A / sigma.inv;
        if (min(B) > r[1] && max(B) < r[2]) {
            break;
        }
        if (h==1) {
            lagged.cov[,,1] <- A;
            lagged.cor[,,1] <- B;
        } else {
            lagged.cov <- abind(lagged.cov, A, along=3);
            lagged.cor <- abind(lagged.cor, B, along=3);
        }
        h <- h + 1;
    }
    return (list(lagged.cov=lagged.cov, lagged.cor=lagged.cor));
}

largeComp <- function(data, level) {
    return(data * (abs(data) > level));
}

hillEstimate <- function(X, prob=0.97) {
    q <- quantile(X, prob);
    if (prob > 0.5) {
        return(1/mean(log(X[which(X > q)]/q)));
    } else {
        return(1/mean(log(X[which(X < q)]/q)));
    }
}

hillPlot <- function(X, prob=0.95, view=c(0.5, 10))
{
    T <- sort(X);
    l <- length(T);
    m <- floor(l * (1 - prob));
    V <- (l-m):(l-1);
    W <- rep(NA, length(V));
    for (i in V) {
        W[i-l+m+1] <- 1/mean(log(T[i:l]/T[i]));
    }
    K <- rev(l-V);
    alpha <- rev(W);
    plot(K, alpha, type="l", ylim=view);
    grid();
    return (list(K=rev(l-V), alpha=rev(W)));
}

pikandsEstimate <- function(X, prob=0.97) {
    q <- quantile(X, prob);
    Y <- sort(X, decreasing=TRUE);
    for (k in (1:length(X))) {
        if (Y[k] > q) {
            k <- k + 1;
        } else {
            break;
        }
    }
    return(log(2)/log((Y[k] - Y[2*k])/(Y[2*k] - Y[4*k])));
}

dekkersEinmahldeHaan <- function(X, prob=0.97) {
    q <- quantile(X, prob);
    Y <- log(X[which(X > q)]);
    H1 <- mean(Y - log(q));
    H2 <- mean((Y - log(q))^2);
    result <- 1 + H1 + 1/2/(H1^2/H2 - 1);
    return(1/result);
}

objf.garch <- function(vartheta, eps,n,sig2init,petit=sqrt(.Machine$double.eps),r0=10){
  omega <- vartheta[1]
  alpha <- vartheta[2]
  beta <- vartheta[3]
  sig2<-rep(0,n)
  sig2[1]<-sig2init
  for(t in 2:n){
    sig2[t]<-omega+alpha*eps[t-1]^2+beta*sig2[t-1]
  }
  qml <- mean(eps[(r0+1):n]^2/sig2[(r0+1):n]+log(sig2[(r0+1):n]))
  qml }
#
VarAsymp<- function(omega,alpha,beta,eps,sig2init,petit,r0=10){
  n <- length(eps)
  dersigma2<-matrix(0,nrow=3,ncol=n)
  sig2<-rep(0,n)
  sig2[1]<-sig2init
  for(t in 2:n){
    vec<-c(1,eps[t-1]^2,sig2[t-1])
    sig2[t]<-omega+beta*sig2[t-1]+alpha*eps[t-1]^2
    dersigma2[1:3,t]<-vec/sig2[t]+beta*dersigma2[1:3,(t-1)]
  }
  eta <- eps[(r0+1):n]/sqrt(sig2)[(r0+1):n]
  eta <- eta/sd(eta)

  J<-dersigma2[1:3,(r0+1):n]%*%t(dersigma2[1:3,(r0+1):n])/(n-r0)
  kappa4<-mean(eta^4)

{if(kappa(J)<1/petit) inv<-solve(J) else inv<-matrix(0,nrow=3,ncol=3)}
var<-(kappa4-1)*inv
list(var=var,residus=eta)
}


estimGARCH<- function(omega,alpha,beta,eps,petit=sqrt(.Machine$double.eps),r0=10)
{
  valinit<-c(omega,alpha,beta)
  n <- length(eps)
  sig2init<-var(eps[1:min(n,5)])
  res <- nlminb(valinit,objf.garch,lower=c(petit,0,0),
                upper=c(Inf,1,1), eps=eps,n=n,sig2init=sig2init)
  omega <- res$par[1]
  alpha<- res$par[2]
  beta <- res$par[3]
  var<-VarAsymp(omega,alpha,beta,eps,sig2init,petit=sqrt(.Machine$double.eps),r0=10)
  list(coef=c(omega,alpha,beta),residus=var$residus,var=var$var)
}
