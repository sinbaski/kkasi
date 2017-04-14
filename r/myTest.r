rm(list=ls());
source("libxxie.r")
library(parallel)
library(doParallel);
library(foreach);

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

data <- getInterpolatedReturns("2010-01-01", "2015-01-01",
                            tables=tables, suffix="_US");
names <- data$assets;
X <- data$ret;

p <- dim(X)[2];
n <- dim(X)[1];
k <- max(floor(dim(X)[1]*0.03), 50);

alpha <- foreach (i = 1:p, .combine=rbind) %do% {
    hillEstimate(-X[, i], k)
}

pdf("../papers/FX/HillTest_Information_Technology.pdf")
plot(1, 1, type="n", xlim=c(1, p), ylim=c(1, p),
     xaxt="n", yaxt="n", xlab="", ylab="",
     main="Information Technology");
for (i in 1:p) {
    for (j in i:p) {
        sd <- sqrt(alpha[i]^2 + alpha[j]^2) / sqrt(k);
        if (i == j) {
            color = "black";
        } else if (abs(alpha[i] - alpha[j]) <
                   qnorm(mean=0, sd=sd, p=0.85)) {
            color="grey";
        } else if (abs(alpha[i] - alpha[j]) <
                   qnorm(mean=0, sd=sd, p=0.9)) {
            color="#00FF00";
        } else if (abs(alpha[i] - alpha[j]) <
                   qnorm(mean=0, sd=sd, p=0.95)) {
            color="#0000FF";
        } else {
            color="#FF0000";
        }
        points(x=i, y=j, pch=19, cex=2, col=color);
    }
}
axis(side=1, at=1:p, labels=gsub("_US", "", names), las=2);
axis(side=2, at=1:p, labels=gsub("_US", "", names), las=1);
## df <- 2 + 0.1 * (1:length(names));
## axis(side=1, at=1:p, labels=df, las=2);
## axis(side=2, at=1:p, labels=df, las=1);
abline(h=1:p, lty=3, col="gray");
abline(v=1:p, lty=3, col="gray");
dev.off();
