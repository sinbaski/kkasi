library(RMySQL);
rm(list=ls());

n.obs = 1000;
last.day = "2015-05-22";
database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host=Sys.getenv("PB"));

assetSet <- "indices";
results = dbSendQuery(database, sprintf("select tblname from %s;", assetSet));
tables <- fetch(results, n=-1)[[1]];
p = length(tables);
dbClearResult(results);

R <- matrix(nrow=n.obs, ncol=length(tables));
for (j in 1 : length(tables)) {
    results <- dbSendQuery(database, sprintf("select closing from %s
where day <= '%s' order by day desc limit %d;", tables[j], last.day, n.obs+1));
    prices <- rev(log(fetch(results, n=-1)[[1]]));
    dbClearResult(results);
    R[,j] <- tail(prices, -1) - head(prices, -1);
}
dbDisconnect(database);

C <- t(R) %*% R / n.obs;
E <- eigen(C, symmetric=TRUE);
M <- E$vectors;
for (k in 1 : dim(M)[1]) {
    M[,k] <- M[,k]/sum(abs(M[,k]));
}
factors <- R %*% M;
corr <- matrix(nrow=30, ncol=dim(M)[1]);
for (k in 1:dim(M)[1]) {
    A <- acf(factors[,k], plot=FALSE)$acf;
    corr[,k] <- A[-1];
}
lm1 <- lm(R[,1] ~ factors);



    
