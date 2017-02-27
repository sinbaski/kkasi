rm(list=ls());
source("libxxie.r");

## Consumer staples
## tables <- c(
##     "ADM",
##     "BF_series_B",
##     "CAG",
## ##    "CL",
##     "CLX",
##     "COST",
##     "CPB",
##     "CVS",
##     "DPS",
## ##    "EL",
##     "GIS",
##     "HRL",
##     "HSY",
##     "K",
##     "KMB",
##     "KO",
##     "KR",
##     "MDLZ",
##     "MJN",
##     "MKC",
##     "MNST",
##     "MO",
##     "PEP",
##     "PG",
##     "PM",
##     "RAI",
##     "SJM",
##     "STZ",
##     "SYY",
##     "TAP",
##     "TSN",
## ##    "WBA",
## ##    "WFM",
##     "WMT"
## );

#Energy
## tables <- c (
##     "APA",
##     "APC",
##     "BHI",
##     "CHK",
##     "COG",
##     "COP",
##     "DO",
##     "DVN",
##     "EOG",
##     "EQT",
##     "FTI",
##     "HAL",
##     "HES",
##     "HP",
##     "KMI",
##     "MPC",
##     "MRO",
##     "MUR",
##     "NBL",
##     "NFX",
##     "NOV",
##     "OKE",
##     "OXY",
##     "PSX",
##     "PXD",
##     "RIG",
##     "RRC",
##     "SE",
##     "SLB",
##     "SWN",
##     "TSO",
##     "VLO",
##     "WMB",
##     "XEC",
##     "XOM"
## );

## Information Technology
tables <- c(
    "ADBE",
    "ADI",
    "ADP",
    "ADSK",
    "AKAM",
    "AMAT",
    "CA",
    "CSCO",
    "CTSH",
    "CTXS",
    "EA",
    "EBAY",
    "FFIV",
    "FISV",
    "HPQ",
    "HRS",
    "IBM",
    "INTC",
    "INTU",
    "JNPR",
    "KLAC",
    "LLTC",
    "LRCX",
    "MCHP",
    "MSFT",
    "MSI",
    "MU",
    "NTAP",
    "NVDA",
    "ORCL",
    "PAYX",
    "QCOM",
    "RHT",
    "SWKS",
    "SYMC",
    "TSS",
    "TXN",
    "VRSN",
    "WDC",
    "XLNX",
    "XRX",
    "YHOO"
);


data <- getInterpolatedPrices(day1="2000-01-01",
                               day2="2015-01-01",
                               suffix="_US",
                               tables=tables
                               ## assetSet="SP500_components"
                               );
price <- data$P;
assets <- data$assets;

X <- diff(log(price));
n <- dim(price)[1];
for (i in 1:dim(X)[2]) {
    I <- which(price[1:(n-1), i]/price[2:n, i] > 1.5);
    if (length(I) > 0) {
        M <- round(price[I, i]/price[I+1, i]);
        X[I, i] <- log(price[I+1, i]) - log(price[I]/M);
    }
    
    J <- which(price[2:n, i]/price[1:(n-1), i] > 1.5);
    if (length(J) > 0) {
        M <- round(price[J+1, i]/price[J, i]);
        X[J, i] <- log(price[J+1, i]) - log(price[J]*M);
    }
}

x.k <- round(log(dim(X)[1])^2);

tail.x <- matrix(NA, nrow=3, ncol=dim(X)[2]);
A.x <- matrix(NA, nrow=2, ncol=dim(X)[2]);
pickands.x.up <- matrix(NA, ncol=dim(X)[2], nrow=2);
pickands.x.down <- matrix(NA, ncol=dim(X)[2], nrow=2);
for (i in 1:dim(X)[2]) {
    tail.x[1, i] <- hillEstimate(X[, i]);
    tail.x[2, i] <- hillEstimate(-X[, i]);
    tail.x[3, i] <- 2 * hillEstimate(X[, i]^2);    

    A.x[1, i] <- scaleEstimate(X[, i], x.k);
    A.x[2, i] <- scaleEstimate(-X[, i], x.k);

    result <- pickandsEstimate(X[, i]);
    pickands.x.up[1, i] <- result$xi;
    pickands.x.up[2, i] <- result$k;

    result <- pickandsEstimate(-X[, i]);
    pickands.x.down[1, i] <- result$xi;
    pickands.x.down[2, i] <- result$k;
}

## Y <- matrix(rt(df=mean(tail.x),
##                n=length(X)),
##             nrow=dim(X)[1],
##             ncol=dim(X)[2]);

## A.y <- apply(Y, MARGIN=2,
##              FUN=function(S) {
##                  scaleEstimate(S, x.k)
##              });

## tail.y <- apply(Y, MARGIN=2,
##                 FUN=hillEstimate);

## Hoga.x <- apply(X, MARGIN=2,
##                 FUN=function(S) {
##                     max(HogaTest(S, p=3.0e-2, t0=0.1))
##                 });
## Hoga.y <- apply(Y, MARGIN=2,
##                 FUN=function(S) {
##                     max(HogaTest(S, p=3.0e-2, t0=0.1))
##                 });

## Quintos.x <- apply(X, MARGIN=2,
##                 FUN=function(S) {
##                     max(QuintosFanRollingTest(S, k=2.5e-2, gam0=0.1,
##                                               variant="GARCH"))
##                 });
## Quintos.y <- apply(Y, MARGIN=2,
##                 FUN=function(S) {
##                     max(QuintosFanRollingTest(S, k=2.5e-2, gam0=0.1,
##                                               variant="iid"))
##                 });
## Hn <- asymptoticDist(n.paths=2000, n.steps=1000, t0=0.2);
## Qn <- QuintosFanRollingDist(gam0=0.1, n.paths=2000, n.steps=1000);


J <- which(abs(tail.x[1, ] - tail.x[2, ]) < 0.5);
q <- qnorm(0.975);
M <- tail.x[2, ] + tail.x[2, ]/sqrt(x.k) * q;
m <- tail.x[2, ] - tail.x[2, ]/sqrt(x.k) * q;
graphics.off();
## par(mfrow=c(1, 2));
pdf(file="../papers/FX/Information_Technology_lower.pdf",
    width=7, height=3.5);
plot(tail.x[2, ],
     main=expression(alpha[L]),
     ylim=c(min(m), max(M)),
     t="p",xaxt="n", xlab="", ylab="",
     );

## for (i in J) {
for (i in 1:dim(X)[2]) {
    x <- i + c(-0.5, 0.5);
    y <- c(m[i], M[i]);
    co <- expand.grid(x, y);
    polygon(x=co[, 1], y=co[, 2], col="#888888");
}
axis(side=1, at=1:length(tables), labels=gsub("_series_", ".", tables),
     las=2);
abline(v=1:length(tables), lty=3);
abline(h=c(median(m), median(M)), lty=1, col="#FF0000", lwd=2);
dev.off();


q <- qnorm(0.975);
sigma <- sqrt(pickandsVar(pickands.x[1, ]) / pickands.x[2, ]);
M <- pickands.x[1, ] + sigma * q;
m <- pickands.x[1, ] - sigma * q;
graphics.off();
## par(mfrow=c(1, 2));
pdf(file="../papers/FX/Information_Technology_Pickands.pdf",
    width=7, height=3.5);
plot(pickands.x[1, ],
     main=expression(alpha),
     ylim=c(min(m), max(M)),
     t="p",xaxt="n", xlab="", ylab="",
     );

for (i in 1:dim(X)[2]) {
    x <- i + c(-0.5, 0.5);
    y <- c(m[i], M[i]);
    co <- expand.grid(x, y);
    polygon(x=co[, 1], y=co[, 2], col="#888888");
}
axis(side=1, at=1:length(tables), labels=gsub("_series_", ".", tables),
     las=2);
abline(v=1:length(tables), lty=3);
abline(h=c(median(m), median(M)), lty=1, col="#FF0000", lwd=2);
dev.off();




lx <- log10(A.x);
ly <- log10(A.y);
l.m <- min(c(lx, ly));
l.M <- max(c(lx, ly));
pdf("../papers/FX/Information_Technology_Hill_scales.pdf")
plot(lx, pch=0, ylim=c(-9, 0), xaxt="n",
     xlab="", main=expression(log[10](A)),
     ylab="");
axis(side=1, at=1:length(tables),
     labels=gsub("_series_", ".", tables),
     las=2);
points(ly, pch=16, col="#FF0000");
abline(v=1:length(tables), lty=3);
abline(h=seq(from=-9, by=0.5, to=0), lty=3);
dev.off();

R <- apply(X,
           MARGIN=2,
           FUN=function(S) {
               result <- myPotEstimate(S, ceiling(length(S) * 0.03));
               return(result);
           });
params <- matrix(NA, nrow=dim(X)[2], ncol=3);
for (i in 1:dim(X)[2])
{
    params[i, 1] <- R[[i]]$index;
    params[i, 2] <- R[[i]]$scale;
    params[i, 3] <- R[[i]]$shift;
}

## pdf("../papers/FX/Information_Technology_OLS_estimates.pdf");
## par(mfrow=c(3, 1));
## plot(params[-21, 1], xlab="", ylab=expression(alpha), xaxt="n");
## axis(side=1, at=1:length(tables),
##      labels=gsub("_series_", ".", tables),
##      las=2);
## abline(v=1:(dim(X)[2]-1), lty=3);
## abline(h=seq(from=1.0, by=0.5, to=9), lty=3);

## plot(log10(params[-21, 2]), xlab="", ylab=expression(log[10](A)), xaxt="n")
## axis(side=1, at=1:length(tables),
##      labels=gsub("_series_", ".", tables),
##      las=2);
## abline(v=1:(dim(X)[2]-1), lty=3);
## abline(h=seq(from=-10, by=0.5, to=-5), lty=3);

## plot(params[-21, 3], xlab="", ylab="shift", xaxt="n");
## axis(side=1, at=1:length(tables),
##      labels=gsub("_series_", ".", tables),
##      las=2);
## abline(v=1:(dim(X)[2]-1), lty=3);
## dev.off();

pdf("../papers/FX/Information_Technology_OLS_estimates.pdf");
par(mfrow=c(3, 1));
plot(params[, 1], xlab="", ylab=expression(alpha), xaxt="n");
axis(side=1, at=1:length(tables),
     labels=gsub("_series_", ".", tables),
     las=2);
abline(v=1:(dim(X)[2]), lty=3);
abline(h=seq(from=1.0, by=0.5, to=10), lty=3);

plot(log10(params[, 2]), xlab="", ylab=expression(log[10](A)), xaxt="n")
axis(side=1, at=1:length(tables),
     labels=gsub("_series_", ".", tables),
     las=2);
abline(v=1:(dim(X)[2]), lty=3);
abline(h=seq(from=-15, by=0.5, to=-4), lty=3);

plot(params[, 3], xlab="", ylab="shift", xaxt="n");
axis(side=1, at=1:length(tables),
     labels=gsub("_series_", ".", tables),
     las=2);
abline(v=1:(dim(X)[2]), lty=3);
dev.off();


## Comparison of upper- and lower-tail indices
plot(tail.x[1, ], t="b", col="#00FF00", ylim=c(min(tail.x), max(tail.x)));

lines(tail.x[2, ], col="#FF0000");
points(tail.x[2, ], col="#FF0000");

lines(tail.x[3, ], col="#000000");
points(tail.x[3, ], col="#000000");

abline(h=seq(from=2.5, to=4.5, by=0.5), lty=3);
abline(v=1:dim(X)[2], lty=3);
