rm(list=ls());
source("libxxie.r")
library(parallel)
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

## ## Consumer staples
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
names <- data$assets;
X <- data$ret;

Fn <- asymptoticDist(1000, 2000, t0=0.1);

## A <- matrix(NA, nrow=100, ncol=dim(X)[2]);
## for (i in 1:length(names)) {
##     A[, i] <- HogaTest(-X[, i], p=0.02, 0.1);
## }

## Simulate from t-distributions
## X <- matrix(NA, nrow=dim(data$ret)[1], ncol=dim(data$ret)[2]);
## for (i in 1:dim(X)[2]) {
##     X[, i] <- rt(n=dim(X)[1], df=2+0.1*i);
## }

B <- matrix(NA, nrow=dim(X)[2], ncol=dim(X)[2]);
for (i in 1:(length(names)-1)) {
    myfun <- function(j) {
        h <- NA;
        tryCatch( {
            h <- max(HogaTest(-c(X[, c(i, j)]), p=0.02, 0.1))
        }, error=function(e) h <- NA
        );
        return(h);
    }
    ## B[i, (i+1):dim(X)[2]] <- unlist(lapply((i+1):length(names), myfun));
    B[i, (i+1):dim(X)[2]] <- unlist(mclapply((i+1):length(names), myfun, mc.cores=detectCores()));
    print(paste("Row ", i));
}


p <- dim(X)[2];
pdf("../papers/FX/Hoga_Energy_pair.pdf")
## pdf("../papers/FX/t_sim_pair.pdf")
plot(1, 1, type="n", xlim=c(1, p), ylim=c(1, p),
     xlab="", ylab="",
     xaxt="n", yaxt="n", main="Energy");
for (i in 1:p) {
    for (j in i:p) {
        if (is.na(B[i, j])) {
            color = "black";
        } else if (B[i, j] == 0) {
            color="white";
        } else if (B[i, j] < quantile(Fn, 0.85)) {
            color="grey";
        } else if (B[i, j] < quantile(Fn, 0.9)) {
            color="#00FF00";
        } else if (B[i, j] < quantile(Fn, 0.95)) {
            color="#0000FF";
        } else {
            color="#FF0000";
        }
        points(x=i, y=j, pch=19, cex=2, col=color);
    }
}

labels <- gsub("_series_", ".", gsub("_US", "", names));
axis(side=1, at=1:p, labels=labels, las=2);
axis(side=2, at=1:p, labels=labels, las=1);
## df <- 2 + 0.1 * (1:length(names));
## axis(side=1, at=1:p, labels=df, las=2);
## axis(side=2, at=1:p, labels=df, las=1);
abline(h=1:p, lty=3, col="gray");
abline(v=1:p, lty=3, col="gray");
dev.off();

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


## p <- 1000;
## n <- dim(X)[1];
## t0 <- 0.1;
## T <- matrix(rt(n=n*p, df=3), nrow=n, ncol=p);
## H <- rep(NA, p);
## for (i in 281:300) {
##     H[i] <- max(HogaTest(T[, i], 0.02, t0));
##     print(paste("series ", i, ": ", H[i]));
## }
## write.table(H[281:300], file="../papers/FX/T3_Hoga_Statistics.dat",
##             quote=F, col.names=F, row.names=F,
##             append=TRUE
##             );
F <- read.table(file="../papers/FX/HogaAsymptotic.dat")$V1;
H <- read.table(file="../papers/FX/T3_Hoga_Statistics.dat")$V1;
G <- read.table(file="../papers/FX/T4_Hoga_Statistics.dat")$V1;

f <- density(F);
h <- density(H);
g <- density(G);

pdf("../papers/FX/Hoga_AsymptoticDistribution.pdf");
plot(
    f$x, f$y, type="l",
    xlab=expression(Q[1]),
    ylab="density",
    xlim=c(0, 200),
    ylim=c(0, max(c(f$y, h$y, g$y))),
    lwd=2
    );
lines(h$x, h$y, col="red", lwd=1.5);
lines(g$x, g$y, col="blue", lwd=1.5);
legend(
    "topright",
    lwd=c(2, 1.5, 1.5), col=c("black", "red", "blue"),
    legend=c(
        "Asymptotic",
        expression(t(3)),
        expression(t(4)))
);
dev.off();

## R <- mclapply(
##     X=1:p,
##     mc.cores=detectCores(),
##     FUN=function(i) {
##         return(max(HogaTest(T[, i], 0.02, t0)));
##     }
## );
## R <- unlist(R);
## write.table(R, file="../papers/FX/T3_Hoga_Statistics.dat",
##             quote=F, col.names=F, row.names=F,
##             append=TRUE

## H <- apply(
##     T, MARGIN=2, FUN=function(V) {
##         return(Hogatest(V, 0.02, t0));
##     }
## );



