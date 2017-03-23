rm(list=ls());
require(MASS);
require(gplots);
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

X <- getInterpolatedReturns("2010-01-01", "2015-01-01",
                            tables=tables, suffix="_US");
X <- X$ret;
my.den <- function(x, K, alpha)
{
    return(alpha * K^alpha / (K - x)^(alpha + 1));
}
params <- matrix(NA, ncol=2, nrow=dim(X)[2]);
for (i in 1:dim(X)[2]) {
    R <- X[, i];
    hill <- hillEstimate(-R);
    if (is.na(hill)) next;
    tryCatch( {
        result <- fitdistr(R[R < 0],
                           function(x, K) my.den(x, K, hill),
                           list(K=1),
                           lower=c(1.0e-8),
                           upper=c(Inf));
        params[i, 1] <- result$estimate;
        params[i, 2] <- hill;
    }, error = function(e) e);
}

I <- !is.na(params[, 1]);
alpha <- params[I, 2];
K <- params[I, 1];

sd = K * sqrt(1 + 2 / alpha);

f <- function(X) {
    X[X < 0] <- 0;
    return(X);
}

pdf(file="../papers/FX/Information_Technology_alpha_K_ci.pdf");
A <- sort(alpha, index.return=TRUE);
plot(alpha[A$ix], (K + 2*sd)[A$ix],
     type="l", ylim=c(0, max(K + 2*sd)),
     xlab=expression(alpha),
     ylab=expression(K),
     main="Information Technology");
polygon(
    x=c(alpha[A$ix], rev(alpha[A$ix])),
    y=c(K[A$ix] + 2*sd[A$ix], rep(0, length(K))),
    col="grey"
);
plotCI(alpha, K,
     uiw=0, liw=0,
     ui=K + 2*sd,
     li=rep(0, length(K)),
     col="#FF0000",
     lwd=1.5,
     xlab=expression(alpha),
     ylab=expression(K),
     add=T,
     main="Information Technology");
grid(nx=20);
dev.off();



