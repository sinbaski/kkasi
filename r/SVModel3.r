rm(list=ls());
library(RMySQL);
source("libxxie.r");

database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host="localhost");

sector="Energy";

tables <- c(
	"APA_US",
	"APC_US",
	"BHI_US",
	"CAM_US",
	"CHK_US",
	"CNX_US",
	"COG_US",
	"COP_US",
	"CVX_US",
	"DO_US",
	"DVN_US",
	"EOG_US",
	"EQT_US",
	"ESV_US",
	"FTI_US",
	"HAL_US",
	"HES_US",
	"HP_US",
	"KMI_US",
	"MPC_US",
	"MRO_US",
	"MUR_US",
	"NBL_US",
	"NFX_US",
	"NOV_US",
	"OKE_US",
	"OXY_US",
	"PXD_US",
	"RIG_US",
	"RRC_US",
	"SE_US",
	"SLB_US",
	"SWN_US",
	"TSO_US",
	"VLO_US",
	"WMB_US",
	"XEC_US",
	"XOM_US"
);
ret <- getAssetReturns("2010-04-28", "2015-08-06", tables, "localhost");
p <- dim(ret)[2];
ratios <- rep(NA, 5);
for (k in 1:length(ratios))
{
    C <- cov(ret[1:(k*200),]);    
    E <- eigen(C, only.values=TRUE);
    ratios[k] <- (E$values[1] - E$values[p])/sum(diag(C));
}
pdf(sprintf("~/Documents/%s.pdf", sector));
plot((1:5)*200, ratios, type="b", xlab="n", 
     ylab=expression((lambda[1]-lambda[p])/"trace"),
     main=sector);
dev.off();
dbDisconnect(database);


