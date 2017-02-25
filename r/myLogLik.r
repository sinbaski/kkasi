rm(list=ls());
require(MASS);
source("libxxie.r")
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
X <- getInterpolatedReturns("2000-01-01", "2015-01-01",
                            tables=tables, suffix="_US");
X <- X$ret;
my.den <- function(x, K, alpha)
{
    return(alpha * K^alpha / (K - x)^(alpha + 1));
}
params <- matrix(NA, ncol=2, nrow=dim(X)[2]);
for (i in 1:dim(X)[2]) {
    R <- X[, i];
    R <- R[R < 0];
    hill <- hillEstimate(-R);
    tryCatch( {
        result <- fitdistr(R[R < 0], my.den,
                           list(K=1, alpha=hill),
                           lower=c(1.0e-8, 1),
                           upper=c(Inf, Inf));
        params[i, ] <- result$estimate;
    }, error = function(e) e);
}



