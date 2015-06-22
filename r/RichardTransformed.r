rm(list=ls());
library(RMySQL);
source("../../cake/libxxie.r");

day2 = '2015-02-28';
day1 = '2010-01-01';

assetSet <- "SP500_components";

database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host=Sys.getenv("PB"));
results = dbSendQuery(database, sprintf("select symbol from %s;",
    assetSet));
tables <- fetch(results, n=-1)[[1]];
dbClearResult(results);
p <- length(tables);

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
prices <- matrix(NA, nrow=n.records, ncol=sum(to.include));

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
R.trans <- matrix(NA, nrow=dim(R)[1]-1, ncol=dim(R)[2]);

for (i in 1 : length(stocks.included)) {
    ## Fn <- ecdf(R[, i]);
    ## U <- Fn(R[, i]);
    ## U <- U[U < 1];
    R.trans[, i] <- -1/log(head(rank(R[,i]), -1)/n.records);
}
p <- sum(to.include);
E <- eigen((n.records * p)^(-2) * t(R.trans) %*% R.trans);
lambda <- sort(E$values);

plot(lambda[1:p-1]/lambda[2:p]);


