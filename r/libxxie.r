## library(abind);
library(parallel)
## library(rmgarch);

###
# A matrix of return values derived from price data on common days
# of all stocks in the given tables.
### 
getAssetReturns <- function(day1, day2, tables, lag,
                            col.name, host)
{
    R <- getAssetPrices(day1, day2, tables, lag, col.name, host);
    return(diff(log(R)));
}

getAssetPrices <- function(day1, day2, tables, lag, col.name, host)
{
    require(RMySQL);
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
        R[,i] <- prices[I];
    }
    dbDisconnect(database);
    return (R);
}

computeDeviations <- function(X, period)
{
    n <- length(X);
    k <- ceiling(n/period);
    dev <- rep(NA, n);
    for(i in period:n) {
        days <- 1 : period;
        line <- lm(X[(i - period + 1) : i] ~ days);
        dev[i] <- line$residuals[period];
    }
    dev <- dev[period:n];
    
    ## for (i in 1:k) {
    ##     offset <- (i - 1) * period;
    ##     l <- min(period, n - offset);
    ##     days <- 1:l;
    ##     line <- lm(X[(offset + 1) : (offset + l)] ~ days);
    ##     dev[(offset + 1) : (offset + l)] <- line$residuals;
    ## }
    return(dev);
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


getInterpolatedPrices <- function(day1, day2, assetSet="",
                                   tables="", suffix="")
{
    require(RMySQL);
    database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
        dbname='avanza', host="localhost");
    if (nchar(assetSet) > 0) {
        results = dbSendQuery(
            database,
            sprintf("select symbol from %s;", assetSet)
        );
        tables <- fetch(results, n=-1)[[1]];
        dbClearResult(results);
        suffix <- "";
        if (assetSet == "SP500_components") {
            suffix <- "_US";
        }
    }
    n.stocks <- length(tables);

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
    return(list(P=prices, assets=tables[to.include]));
}

getInterpolatedReturns <- function(day1, day2, assetSet="",
                                   tables="", suffix="")
{
    R <- getInterpolatedPrices(day1, day2, assetSet,
                               tables, suffix);
    return(list(ret=diff(log(R$P)), assets=R$assets));
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

hillEstimate <- function(X, k=round(log(length(X))^2)) {
    X.sorted <- sort(X, decreasing=TRUE);
    return(1/mean(log(X.sorted[1:k]/X.sorted[k+1])));
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
    plot(K, alpha, type="l", ylim=view,
         ylab=expression(alpha),
         xlab="Upper Order Statistic"
         );
    grid();
    return (list(K=rev(l-V), alpha=rev(W)));
}

## pickandsEstimate <- function(X, prob=0.97) {
##     q <- quantile(X, prob);
##     Y <- sort(X, decreasing=TRUE);
##     for (k in (1:length(X))) {
##         if (Y[k] > q) {
##             k <- k + 1;
##         } else {
##             break;
##         }
##     }
##     return(log(2)/log((Y[k] - Y[2*k])/(Y[2*k] - Y[4*k])));
## }

## The procedure programmed here is only meant to
## estimate the extreme value index when it is postiive,
## i.e. when the distribution is in the MDA of Frechet.
pickandsEstimate <- function(X, k=NA)
{
    X.sorted <- sort(X, decreasing=TRUE);
    if (!is.na(k)) {
        a <- X.sorted[k] - X.sorted[2 * k];
        b <- X.sorted[2 * k] - X.sorted[4 * k];
        xi <- log(a/b) / log(2);
        return(xi);
    }
    ## mef <- meanExcess(X);
    ## M <- which.max(mef$mef);
    ## k <- length(mef$thresholds) + 1 - M + 1;
    k <- 1;
    L <- round(length(X) * 0.03);
    xi <- rep(NA, L - k + 1);
    for (i in k:L) {
        a <- X.sorted[i] - X.sorted[2 * i];
        b <- X.sorted[2 * i] - X.sorted[4 * i];
        xi[i-k+1] <- log(a/b) / log(2);
    }
    xi <- sort(xi[xi > 0]);
    k <- floor(length(xi)/2);
    return(list(k=k, xi=xi[k]));
}

pickandsVar <- function(xi)
{
    r <- xi^2 * (2^(2*xi + 1) + 1);
    t <- 2 * log(2) * (2^xi - 1);
    return(r/t^2);
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

weissmanEstimate <- function (X, s, t, p)
{
    n <- length(X);
    k <- max(50, n * 0.03);
    X.sorted <- sort(X[ceiling(n*s):floor(n*t)], decreasing=TRUE);
    N <- ceiling(k * (t - s));
    hill.inv <- mean(log(X.sorted[1:(N-1)] / X.sorted[N]));
    return((n * p / k)^(-hill.inv) * X.sorted[N]);
}

scaleEstimate <- function(X, m)
{
    n <- length(X);
    X.sorted <- sort(X, decreasing=TRUE);
    X.m <- X.sorted[m+1];
    hill <- 1/mean(log(X.sorted[1:m] / X.sorted[m+1]));
    return(X.m^hill * m/n);
}

meanExcess <- function(Z, U=NA)
{
    if (is.na(U)) {
        Z <- sort(Z);
        n <- length(Z);
        k <- n - round(log(length(Z))^2) + 1;
        R <- sapply(k:(n-1), function(i) {
            mean(Z[(i+1):n] - Z[i])
        });
        U <- Z[k:(n-1)];
        return(list(thresholds=U, mef=R));
    }
    i <- 1;
    ME <- rep(NA, length(U));
    for (u in U) {
        X <- Z[Z > u];
        ME[i] <- mean(X - u);
        i <- i + 1;
    }
    return(list(thresholds=U, mef=ME));
}

myPotEstimate <- function(Z, k)
{
    Z <- sort(Z, decreasing=T);
    U <- rev(Z[2:k]);
    mef <- meanExcess(Z, U);
    n2 <- which.max(mef);
    n1 <- which.min(mef[1:n2]);

    X <- U[n1:n2];
    Y <- mef[n1:n2];
    mdl <- lm(Y ~ X);
    C <- coef(mdl);
    alpha <- 1/C[2] + 1;
    a <- C[1] * (alpha - 1);
    
    ## u <- Z[n1];
    ## Y <- sort(Z[Z > u], decreasing=T);
    ## n <- length(Y);
    p <- rev(1:(k-1))/length(Z);
    p <- log(p[n1:n2]);
    q <- log(U[n1:n2] + a);
    model <- lm(p ~ q);
    A <- exp(coef(model)[1]);
    return(list(index=1/C[2] + 1,
                shift=a,
                scale=A)
           );
    ## u <- Z[k];
    ## fun <- function(b) {
    ##     Y <- log(X/b + 1);
    ##     result <- mean(Y + 1) * mean(exp(-Y)) - 1;
    ##     return(result);
    ## }
    ## Maximum likelihood estimate
    ## This often doesn't work
    ## fu <- fun(u);
    ## if (fu < 0) {
    ##     delta <- 1.5;
    ## } else {
    ##     delta <- 1/1.5;
    ## }
    ## b <- u * delta;
    ## while (fun(b) * fu > 0) {
    ##     b <- b * delta;
    ## }
    ## root <- uniroot(f=fun, interval=c(u, b));
    ## b <- root$root;
    ## a <- b - u;
    ## alpha <- 1/(mean(log(X + b)) - log(b));


}

HogaTest <- function(X, p, t0)
{
    times <- seq(from=t0, to=1-t0, length.out=100);
    ## G <- rep(NA, length(times));
    f1 <- function(s, t) {
        r <- rep(NA, length(s));
        y <- weissmanEstimate(X, 0, t, p);
        for (i in 1:length(s)) {
            x <- weissmanEstimate(X, 0, s[i], p)/y;
            r[i] <- (log(x) * s[i])^2;
        }
        return(r);
    }
    f2 <- function(s, t) {
        r <- rep(NA, length(s));
        y <- weissmanEstimate(X, t, 1, p);
        for (i in 1:length(s)) {
            x <- weissmanEstimate(X, s[i], 1, p)/y;
            r[i] <- (log(x) * (1 - s[i]))^2;
        }
        return(r);
    }
    T1 <- rep(NA, length(times));
    T2 <- rep(NA, length(times));
    T3 <- rep(NA, length(times));
    G <- rep(NA, length(times));
    for (i in 1:length(times)) {
        T1[i] <- weissmanEstimate(X, 0, times[i], p)/
            weissmanEstimate(X, times[i], 1, p);
        T1[i] <- (times[i] * (1 - times[i])*log(T1[i]))^2;
        if (T1[i] == 0) {
            G[i] <- 0;
            next
        }
        
        T2[i] <- integrate(f=function(s) {f1(s, times[i])},
                        lower=t0, upper=times[i],
                        subdivisions=200L)$value;
        T3[i] <- integrate(f=function(s) {f2(s, times[i])},
                        lower=times[i], upper=1 - t0,
                        subdivisions=200L)$value;
        G[i] <- T1[i] / (T2[i] + T3[i]);
    }
    return(G);
}

## HogaTest <- function(X, t0, k, p)
## {  
##     n <- length(X)                                      
##     AnzAntEx <- length(k)                               
##     k <- floor(k * n)                                   
    
    
##     T_0t  <- array(NA, dim=c(n, AnzAntEx))      
##     Tq_0t <- array(NA, dim=c(n, AnzAntEx))      

##     Y <- sort.int(X[1 : floor(n * t0)],
##                   decreasing = TRUE, method = "quick")
    
##     for(i in (floor(n * t0) + 1) : n){
##         Y <- sort.int(c(Y, X[i]),
##                       decreasing = TRUE, method = "quick")      

##         for(j in (1 : AnzAntEx)){      
            
##             T_0t[i, j] <- 1/(k[j]*(i/n)) *
##                 sum(log(Y[1 : (floor(k[j]*(i/n))+1)] /
##                         Y[floor(k[j]*(i/n))+1])) 
##             Tq_0t[i, j] <- Y[floor(k[j]*(i/n))+1] *
##                 ((n*p) / k[j])^(-T_0t[i, j])                            
##         }
##     } 
    
##     T_t1  <- array(NA, dim=c(n, AnzAntEx))      
##     Tq_t1 <- array(NA, dim=c(n, AnzAntEx))      
    
##     Y    <- sort.int(X[(floor(n * (1 - t0)) + 1) : n],
##                      decreasing = TRUE, method = "quick")
    
##     for(i in floor(n * (1 - t0)) : 1){      
##         Y <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")  
##         for(j in (1 : AnzAntEx)){           
##             T_t1[i, j] <- 1/(k[j]*(1 - (i-1)/n)) *
##                 sum(log(Y[1 : (floor(k[j]*(1 - (i-1)/n))+1)] /
##                         Y[floor(k[j]*(1 - (i-1)/n)) + 1]))
##             Tq_t1[i, j] <- Y[floor(k[j]*(1 - (i-1)/n))+1] *
##                 ((n*p) / k[j])^(-T_t1[i, j])
##         }
##     }

##     ## Calculate SN test statistics
##     ## Save realised suprema for different no. of k's
##     maxQ_T  <- array(NA, dim = c(AnzAntEx))
##     maxQ_Tq <- array(NA, dim = c(AnzAntEx)) 
##     Q_T     <- array(NA, dim = c(n, AnzAntEx))    
##     Q_Tq    <- array(NA, dim = c(n, AnzAntEx))    
##     ## Calculating the self-normalizer "AnzAntEx" versch.
##     ## k's für die "5" Schätzer
##     SN      <- array(NA, dim = c(n, AnzAntEx))
##     ## Calculating the self-normalizer "AnzAntEx" versch.
##     ## k's für die "5" Schätzer
##     SN_T    <- array(NA, dim = c(n, AnzAntEx))
##     l_sq    <- log(k / (n*p))^2

##     ## loop over time      
##     for(i in (floor(n*t0)+1) : floor(n*(1-t0))){
##         ## discretize time "t"
##         t <- i/n
##         Q_Tq[i, ] <- k / l_sq *
##             (t * (1-t) * log(Tq_0t[i, ] / Tq_t1[i, ]))^2 
##         Q_T[i, ]  <- k  * (t * (1-t) * (T_0t[i, ] - T_t1[i, ]))^2
        
##         for(m in 1:AnzAntEx){
##             SN[i, m] <- 1/n * k[m] / l_sq[m] *
##                 sum((((floor(n*t0)+1) : i)/n *
##                      log(Tq_0t[(floor(n*t0)+1) : i, m] / Tq_0t[i, m]))^2)
##             SN[i, m] <- SN[i, m] + 1/n * k[m] / l_sq[m] *
##                 sum((((n-(i : floor(n*(1-t0))) + 1) / n) *
##                      log(Tq_t1[i : floor(n*(1-t0)), m] / Tq_t1[i, m]))^2)
            
##             SN_T[i, m] <- 1/n * k[m] *
##                 sum((((floor(n*t0)+1) : i)/n *
##                      (T_0t[(floor(n*t0)+1) : i, m] - T_0t[i, m]))^2)
##             SN_T[i, m] <- SN_T[i, m] + 1/n * k[m] *
##                 sum((((n-(i : floor(n*(1-t0))) + 1) / n) *
##                      (T_t1[i : floor(n*(1-t0)), m] - T_t1[i, m]))^2)
##         }
##     }
##     ## elementwise division
##     maxQ_T  <- colMaxs( Q_T[(floor(n * t0) + 1) : floor(n * (1-t0)), ] /
##                         SN_T[(floor(n * t0) + 1) : floor(n * (1-t0)), ])
##     ## elementwise division
##     maxQ_Tq <- colMaxs(Q_Tq[(floor(n * t0) + 1) : floor(n * (1-t0)), ] /
##                        SN[(floor(n * t0) + 1) : floor(n * (1-t0)), ])   
##     return(list(TI = maxQ_T, EQ = maxQ_Tq))
## }

quantilsup.SN1 <- function(AnzSim, AnzGitter, t0, Quantile){   # AnzSim = Anzahl Wiener-Pfade, AnzGitter = Diskretisierungspunkte der Wiener-Pfade, Quantile = Vektor der zu berechnenden Quantile
      RealZVseqc  <- rep.int(0, AnzSim)                    # Hier werden die "AnzSim"-vielen Realisationen der Grenz-Zufallsvariablen sup_{t\in[t0,1-t0]}(B(t)-tB(1))^{2} gespeichert
      GitterZVseqc<- rep.int(0, AnzGitter)                 # An Gitterpunkten (in t) wird (B(t)-tB(1))^{2} ausgewertet
      numerator   <- rep.int(0, AnzGitter)
      denom1      <- rep.int(0, AnzGitter)
      denom2      <- rep.int(0, AnzGitter)
      for(j in 1:AnzSim){
            B <- rwiener(end = 1, frequency = AnzGitter)                # Simulation eines Standard-Wiener-Prozesses
            for(i in ceiling(AnzGitter*t0) : floor(AnzGitter*(1-t0))){
                  t <- i/AnzGitter
                  numerator[i] <- (B[i] - t * B[AnzGitter])^2
                  denom1[i]    <- 1/AnzGitter * sum((B[ceiling(AnzGitter*t0) : i] - (ceiling(AnzGitter*t0) : i) / (AnzGitter*t) * B[i])^2)
                  denom2[i]    <- 1/AnzGitter * sum((B[AnzGitter] - B[i : floor(AnzGitter*(1-t0))] - (1 - (i : floor(AnzGitter*(1-t0))) / AnzGitter) / (1-t) * (B[AnzGitter] - B[i]))^2)
                  GitterZVseqc[i] <- numerator[i] / (denom1[i] + denom2[i])
            }
            RealZVseqc[j]<- max(GitterZVseqc)      # save realisation
      }
      return(quantile(RealZVseqc, probs = Quantile))
}     # quantile-function gives required quantiles
    
asymptoticDist <- function(n.paths, n.steps, t0)
{
    W <- matrix(rnorm(n.steps * n.paths),
                nrow=n.steps, ncol=n.paths);
    ## W <- unlist(mclapply(1:n.steps, cumsum, mc.cores=detectCores()));
    W <- apply(W, MARGIN=2, FUN=cumsum);
    f <- function(index) {
        X <- W[, index];
        n <- length(X);
        ds <- 1/n;
        a <- ceiling(t0 * n);
        b <- floor((1 - t0) * n);
        results <- rep(NA, b-a+1);
        for (i in a:b) {
            t <- i/n;
            S <- (a : i)/n;
            I1 <- sum((X[a:i] - (S/t) * X[i])^2)/n;
            S <- (i : b)/n;
            T <- ((-X[i:b] + X[n]) - ((1 - S)/(1 - t))*(X[n] - X[i]))^2;
            I2 <- sum(T)/n;
            denom <- (X[i] - t * X[1])^2;
            results[i-a+1] <- denom/(I1 + I2);
        }
        return(max(results));
    }
    ## samples <- apply(W, MARGIN=2, FUN=f);
    samples <- unlist(mclapply(X=1:n.paths, FUN=f, mc.cores=detectCores()));
    ## samples <- lapply(X=n.paths, FUN=f);
    return(samples);
}

## gam0: ratio of sub-sample size over full sample size
## k: ratio of order statistics over sub-sample size
## h: length(X) * [h, 1 - h] is the interval to test for change of alpha.
## variant: iid or GARCH
QuintosFanRollingTest <- function(X, k, h, variant)
{
    if (variant == "GARCH") {
        X <- X^2;
    }
    T <- length(X);
    t0 <- floor(T * h);
    wt <- t0 - 1;
    gam0 <- h;
    ## wt <- floor(gam0 * T);
    mt <- floor(k * gam0 * T);
    f1 <- gam0 * sqrt(k * T);
    V <- rep(NA, T - 2 * t0);
    hill <- rep(NA, length(V));
    alpha.T <- hillEstimate(X, round(k*T));
    for (i in 1:length(V)) {
        t <- t0 + i - 1;
        Y <- X[(t - wt + 1):t];
        hill[i] <- hillEstimate(Y, round(k*wt));
        V[i] <- (hill[i]/alpha.T - 1) * f1;
        if (variant == "iid") {
            eta <- 1;
        } else if (variant == "GARCH") {
            Y.s <- sort(Y, decreasing=TRUE);
            c.wt <- log(Y / Y.s[mt]);
            c.wt[c.wt < 0] = 0;
            d.wt <- Y > Y.s[mt];
            
            chi <- (2 * hill[i])^2 * (1/mt) *
                sum(c.wt[1:(wt - 1)] * c.wt[2:wt]);
            psi <- sum(c.wt[1:(wt - 1)] * d.wt[2:wt]
                       + c.wt[2:wt] * d.wt[1:(wt- 1)]);
            psi <- psi / mt * hill[i];
            omega <- 2 * sum(d.wt[1:(wt - 1)] * d.wt[2:wt]) / mt;
            eta <- 1 + chi + omega - 2 * psi;
        }
        V[i] <- V[i] / sqrt(eta);
    }
    return(V^2);
}

## h: length(X) * [h, 1 - h] is the interval to test for change of alpha.
QuintosFanRollingDist <- function(h, n.paths, n.steps)
{
    ## samples <- matrix(rnorm(n=n.steps * n.paths),
    ##                   nrow=n.steps, ncol=n.paths);
    ## wt <- floor(n.steps * gam0);
    ## f <- function(X)
    ## {
    ##     W <- c(0, cumsum(X)) / sqrt(n.steps);
    ##     stat <- rep(NA, n.steps - 2 * t0);
    ##     for (i in 1:length(stat)) {
    ##         t <- i + t0;
    ##         stat[i] <- W[t] - W[t - wt] - gam0 * W[length(W)];
    ##         stat[i] <- stat[i]^2;
    ##     }
    ##     return(max(stat));
    ## }
    ## Fn <- ecdf(apply(samples, MARGIN=2, FUN=f));
    ## return(Fn);

    V <- rep(NA, n.paths);
    t0 <- floor(n.steps * h);
    wt <- t0 - 1;
    gam0 <- h;
    for (i in 1 : n.paths) {
        W <- rwiener(end=1, frequency=n.steps);
        I <- t0:(n.steps - t0);
        w.bar <- W[I] - W[I - wt] - gam0 * W[n.steps];
        V[i] <- max(w.bar^2);
    }
    return(ecdf(V));
}

dmydist <- function(q, tail.index)
{
    return(tail.index * q^(-tail.index - 1));
}

pmydist <- function(q, tail.index)
{
    prob <- -q^(-tail.index) + 1;
    return(prob);
}

qmydist <- function(p, tail.index)
{
    q <- (1 - p)^(-1/tail.index);
    return(q);
}

rmydist <- function(n, tail.index)
{
    X <- runif(n);
    X <- qmydist(X, tail.index);
    return(X);
}

colMaxs <- function(X)
{
    return(apply(X, MARGIN=2, FUN=max));
}

rwiener <- function(end, frequency)
{
    W <- cumsum(rnorm(n=frequency, mean=0, sd=sqrt(end/frequency)));
    return(W);
}
