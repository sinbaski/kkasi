rm(list=ls());
require(MASS);
require(gplots);
source("libxxie.r")
#Energy
tables <- c (
    "APA",
    "APC",
    "BHI",
    "CHK",
    "COG",
    "COP",
    "DO",
    "DVN",
    "EOG",
    "EQT",
    "FTI",
    "HAL",
    "HES",
    "HP",
    "KMI",
    "MPC",
    "MRO",
    "MUR",
    "NBL",
    "NFX",
    "NOV",
    "OKE",
    "OXY",
    "PSX",
    "PXD",
    "RIG",
    "RRC",
    "SE",
    "SLB",
    "SWN",
    "TSO",
    "VLO",
    "WMB",
    "XEC",
    "XOM"
);

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

data <- getInterpolatedReturns("2010-01-01", "2015-01-01",
                            tables=tables, suffix="_US");
X <- data$ret;
my.den <- function(x, K, alpha)
{
    return(alpha * K^alpha / (K - x)^(alpha + 1));
}
params <- matrix(NA, ncol=3, nrow=dim(X)[2]);
A <- rep(NA, dim(X)[2]);

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
        params[i, ] <- c(result$estimate, hill, sum(R < 0));
    }, error = function(e) e);
}

I <- !is.na(params[, 1]);
alpha <- params[I, 2];
K <- params[I, 1];
N <- params[I, 3];
p <- N/dim(X)[1];

A <- apply(X, MARGIN=2, FUN=function(V) scaleEstimate(-V, max(floor(length(V)[1]*0.03), 60)));
B <- p*K^alpha;

pdf(file="../papers/FX/Information_Technology_scale.pdf")
plot(1:sum(I), log10(A[I]),
     col="red",
     ylim=c(min(c(log10(A[I]), log10(B))),
            max(c(log10(A[I]), log10(B)))),
     main="Information Technology Scale Estimates",
     xlab="Equity", ylab="Scale",
     xaxt="n",
     pch=16);
axis(side=1, at=1:sum(I),
     labels=gsub("_series_", ".", gsub("_US", "", data$asset[I])),
     las=2);
points(1:sum(I), log10(p*K^alpha), pch=2);
dev.off();

sd = K * sqrt(1 + 2 / alpha) / sqrt(N);
p.sd <- sqrt(p - p^2)/sqrt(dim(X)[1]);

pdf(file="../papers/FX/Energy_p.pdf");
plotCI(1:sum(I), p,
       barcol="#000000",
       col="#FF0000",
       lwd=1,
       ui=p + 2*p.sd,
       li=p - 2*p.sd,
       type="p", pch=16,
       ylim=c(min(p - 2*p.sd), max(p + 2*p.sd)),
       xlab="",
       ylab=expression(P(X < 0)),
       xaxt="n",
       main="Energy");
axis(side=1, at=1:sum(I),
     labels=gsub("_series_", ".", gsub("_US", "", data$asset[I])),
     las=2);
dev.off();

## abline(v=1:sum(I), lty=3, col="grey");
## abline(h=seq(from=min(p), to=max(p), length.out=10), lty=3, col="grey");


pdf(file="../papers/FX/Energy_alpha_K_ci.pdf");
A <- sort(alpha, index.return=TRUE);
plot(alpha[A$ix], K[A$ix],
     type="p", pch=16,
     xlim=c(min(alpha), max(alpha)),
     ylim=c(min(K - 2 * sd), max(K + 2 * sd)),
     col="#FF0000",
     xlab=expression(alpha),
     ylab=expression(K),
     main="Energy");
## polygon(
##     x=c(alpha[A$ix], rev(alpha[A$ix])),
##     y=c(K[A$ix] + sd[A$ix], rep(0, length(K))),
##     col="grey"
## );
plotCI(alpha, K,
     uiw=0, liw=0,
     ui=K + 2*sd,
     li=K - 2*sd,
     barcol="#000000",
     col="red",
     lwd=1,
     xlab=expression(alpha),
     ylab=expression(K),
     add=T);
abline(v=seq(from=floor(min(alpha)), by=0.5, to=ceiling(max(alpha))),
       lty=3, col="grey");
abline(h=seq(from=0, by=0.01, to=ceiling(max(K+2*sd))),
       lty=3, col="grey");
dev.off();

