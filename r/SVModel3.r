rm(list=ls());
library(RMySQL);
source("libxxie.r");

database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host="localhost");

sector="Materials";

tables <- c(
    "AA_US",
    "APD_US",
    "ARG_US",
    "AVY_US",
    "BLL_US",
    "CF_US",
    "DD_US",
    "DOW_US",
    "ECL_US",
    "EMN_US",
    "FCX_US",
    "FMC_US",
    "IFF_US",
    "IP_US",
    "LYB_US",
    "MLM_US",
    "MON_US",
    "MOS_US",
    "NEM_US",
    "NUE_US",
    "OI_US",
    "PPG_US",
    "PX_US",
    "SEE_US",
    "SHW_US",
    "SIAL_US",
    "VMC_US"
);
ret <- getAssetReturns("2010-04-28", "2015-08-06", tables, "localhost");
p <- dim(ret)[2];
ratios <- rep(NA, 7);
for (k in 2:length(ratios))
{
    C <- cov(ret[1:(k*200),]);    
    E <- eigen(C, only.values=TRUE);
    ratios[k] <- (E$values[1] - E$values[p])/sum(diag(C));
}
pdf(sprintf("/tmp/%s.pdf", sector));
plot((1:7)*200, ratios, type="b", xlab="p", ylab="",
     main=expression(frac(lambda[1]-lambda[p], "trace")));
dev.off();
dbDisconnect(database);


