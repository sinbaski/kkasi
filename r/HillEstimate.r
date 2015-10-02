rm(list=ls());
library(RMySQL);

day2 = '2015-02-28';
day1 = '2010-01-01';

assetSet <- "SP500_components";

## database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
##     dbname='avanza', host=Sys.getenv("PB"));
database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host="localhost");
results = dbSendQuery(database, sprintf("select symbol from %s;",
    assetSet));
tables <- fetch(results, n=-1)[[1]];
dbClearResult(results);
n.stocks <- length(tables);

results <- dbSendQuery(
    database,
    sprintf("select count(*) from Calendar where weekday(day) <= 4
and day between '%s' and '%s';", day1, day2));
n.records <- fetch(results, n=1)[1, 1];
dbClearResult(results);

to.include <- rep(TRUE, n.stocks);
for (i in 1 : n.stocks) {
    tables[i] <- gsub("[.]", "_", tables[i]);
    tables[i] <- gsub("-", "_series_", tables[i]);
    tables[i] <- paste(tables[i], "_US", sep="");

    results <- dbSendQuery(
        database,
        sprintf("select count(*) from %s where
day between '%s' and '%s';", tables[i], day1, day2)
    );
    n.traded.days <- fetch(results)[1,1];
    dbClearResult(results);

    results <- dbSendQuery(
        database,
        sprintf("select min(day) from %s;",
                tables[i])
    );
    d <- fetch(results)[1, 1];
    dbClearResult(results);

    if (n.traded.days < 1000 || d > day1) {
        to.include[i] <- FALSE;
        next;
    }

}
p <- sum(to.include);
prices <- matrix(NA, nrow=n.records, ncol=p);

results <- dbSendQuery(
    database,
    sprintf("select day from Calendar where weekday(day) <= 4
and day between '%s' and '%s';", day1, day2)
);
days <- fetch(results, n=-1)[[1]];
dbClearResult(results);

stocks.included <- which(to.include);
for (i in 1:length(stocks.included)) {
    results <- dbSendQuery(
        database,
        sprintf("select day, closing from %s
where day between '%s' and '%s'",
                tables[stocks.included[i]],
                day1, day2)
    );
    A <- fetch(results, n=-1);
    dbClearResult(results);
    I <- days %in% A[[1]];
    prices[which(I), i] <- A[[2]];

    if (is.na(prices[1, i])) {
        ## get the last trading day before the period
        results <- dbSendQuery(
            database,
            sprintf("select max(day) from %s where day <= '%s';",
                    tables[stocks.included[i]], day1)
        );
        d <- fetch(results)[1, 1];
        dbClearResult(results);

        ## get the closing price of the last trading day
        results <- dbSendQuery(
            database,
            sprintf("select closing from %s where day = '%s';",
                    tables[stocks.included[i]], d)
        )
        last.price <- fetch(results)[1, 1];
        dbClearResult(results);

        ## fill in the price up to the first available
        j <- 1;
        while (is.na(prices[j, i])) {
            prices[j, i] <- last.price;
            I[j] <- TRUE;
            j <- j + 1;
        }
    }
    J <- which(!I);
    for (j in J) {
        prices[j, i] <- prices[j-1, i];
    }
}
dbDisconnect(database);
R <- diff(log(prices));
tailIndices = matrix(NA, nrow=dim(R)[2], ncol=2);
for (i in 1:dim(R)[2]) {
    X = R[, i];
    a <- quantile(X, probs=0.03);
    tailIndices[i, 1] <- 1/mean(log(X[which(X < a)]/a));

    b <- quantile(X, probs=0.97);
    tailIndices[i, 2] <- 1/mean(log(X[which(X > b)]/b));
}

## R.trfm <- matrix(NA, nrow=dim(R)[1], ncol=dim(R)[2]);

## for (i in 1 : length(stocks.included)) {
##     ## Fn <- ecdf(R[, i]);
##     ## U <- Fn(R[, i]);
##     ## U <- U[U < 1];
##     R.trfm[, i] <- -1/log(rank(R[,i])/(n.records+1));
## }
## E <- eigen((n.records * p)^(-2) * t(R.trfm) %*% R.trfm);

## plot(lambda[2:p]/lambda[1:p-1], type="b", xlim=c(1, 20), ylim=c(0,1));
a <- min(tailIndices[, 1]);
b <- max(tailIndices[, 1]);
pdf("SP500_tail_indices.pdf")
plot(tailIndices[, 1], tailIndices[, 2], type="p",
     xlab="Lower tail index", ylab="Upper tail index",
     xlim=c(a, b), ylim=c(a, b));
X <- seq(a, b, 0.01);
lines(X, X, col="#FF0000");
dev.off();


