rm(list=ls());
source("libxxie.r")
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

## Consumer staples
tables <- c(
    "ADM",
    "BF_series_B",
    "CAG",
##    "CL",
    "CLX",
    "COST",
    "CPB",
    "CVS",
    "DPS",
##    "EL",
    "GIS",
    "HRL",
    "HSY",
    "K",
    "KMB",
    "KO",
    "KR",
    "MDLZ",
    "MJN",
    "MKC",
    "MNST",
    "MO",
    "PEP",
    "PG",
    "PM",
    "RAI",
    "SJM",
    "STZ",
    "SYY",
    "TAP",
    "TSN",
##    "WBA",
##    "WFM",
    "WMT"
);

## Information Technology
## tables <- c(
##     "ADBE",
##     "ADI",
##     "ADP",
##     "ADSK",
##     "AKAM",
##     "AMAT",
##     "CA",
##     "CSCO",
##     "CTSH",
##     "CTXS",
##     "EA",
##     "EBAY",
##     "FFIV",
##     "FISV",
##     "HPQ",
##     "HRS",
##     "IBM",
##     "INTC",
##     "INTU",
##     "JNPR",
##     "KLAC",
##     "LLTC",
##     "LRCX",
##     "MCHP",
##     "MSFT",
##     "MSI",
##     "MU",
##     "NTAP",
##     "NVDA",
##     "ORCL",
##     "PAYX",
##     "QCOM",
##     "RHT",
##     "SWKS",
##     "SYMC",
##     "TSS",
##     "TXN",
##     "VRSN",
##     "WDC",
##     "XLNX",
##     "XRX",
##     "YHOO"
## );

X <- getInterpolatedReturns("2010-01-01", "2015-01-01",
                            tables=tables, suffix="_US");
names <- X$assets;
X <- X$ret;

Fn <- asymptoticDist(1000, 2000, t0=0.1);

A <- matrix(NA, nrow=100, ncol=dim(X)[2]);
for (i in 1:length(names)) {
    A[, i] <- HogaTest(-X[, i], p=0.02, 0.1);
}

B <- matrix(NA, nrow=dim(X)[2], ncol=dim(X)[2]);
for (i in 1:length(names)) {
    for (j in (i+1):length(names)) {
        B[i, j] <- max(HogaTest(-c(X[, c(i, j)]), p=0.02, 0.1));
        B[j, i] <- B[i, j];
    }
}


## barplot(G, xlim=c(0, ceiling(q)),
##         width=1, space=0.2,
##         horiz=T,
##         main="Test Statistics of Energy",
##         xlab=expression(Q), ylab=""
##         );
## axis(side=2, at=seq(from=0.6, by=1.2, length.out=length(G)),
##      labels=gsub("_US", "", names[-1]),
##      las=1);

## abline(v=q, lwd=2, col="#FF0000");



H <- apply(A, FUN=max, MARGIN=2);
pdf("../papers/FX/Hoga_Consumer_Staples_Single.pdf")
barplot(H, xlim=c(0, max(quantile(Fn, 0.95), max(H))),
        width=1, space=0.2,
        horiz=T,
        main="Test Statistics of Consumer Staples",
        xlab=expression(Q), ylab=""
        );
axis(side=2, at=seq(from=0.6, by=1.2, length.out=length(H)),
     labels=gsub("_US", "", names),
     las=1);
abline(v=quantile(Fn, 0.85), lwd=2, col="#00FF00");
abline(v=quantile(Fn, 0.90), lwd=2, col="#0000FF");
abline(v=quantile(Fn, 0.95), lwd=2, col="#FF0000");

dev.off();



