rm(list=ls());
library(RMySQL);

day2 = '2015-02-28';
day1 = '2010-01-01';

assetSet <- "SP500_components";

## database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
##     dbname='avanza', host=Sys.getenv("PB"));
database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
    dbname='avanza', host="localhost");
results = dbSendQuery(database, sprintf("select t1.symbol, t2.sector from
SP500_components as t1 join company_info as t2 on t1.symbol = t2.symbol;",
    assetSet));
H <- fetch(results, n=-1);
tables <- H[[1]];
sectors <- H[[2]];
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
        sprintf(
            "select day, closing from %s where day between '%s' and '%s'",
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

## E <- eigen((n.records * p)^(-2) * t(R) %*% R);

## plot(lambda[2:p]/lambda[1:p-1], type="b", xlim=c(1, 20), ylim=c(0,1));
a <- min(tailIndices[, 1]);
b <- max(tailIndices[, 1]);

pdf("SP500_tail_indices_colored_by_sector.pdf")
results <- dbSendQuery(database, "select distinct(sector) from company_info;");
sector.names <- fetch(results)[[1]];
dbClearResult(results);
colors <- rainbow(length(sector.names));

## Record the eigenvalues
## lambda <- matrix(NA, nrow=1+length(sector.names), ncol=p);
company.num <- rep(NA, length(sector.names));
for (i in 1:length(sector.names)) {
    I <- which(sectors[to.include] == sector.names[i]);
    ## if (i == 1) {
    ##     plot(tailIndices[I, 1], tailIndices[I, 2],
    ##          col=colors[i], pch=i,
    ##          xlim=c(a, 6.8), ylim=c(a, b),
    ##          xlab="Lower tail index", ylab="Upper tail index");
    ## } else {
    ##     points(tailIndices[I, 1], tailIndices[I, 2], col=colors[i], pch=i);
    ## }

    ## pdf(sprintf("../papers/%s_tail_indices.pdf", sector.names[i]));
    ## plot(tailIndices[I, 1], tailIndices[I, 2],
    ##      col="#0000FF", pch=i,
    ##      xlim=c(a, b), ylim=c(a, b),
    ##      main=sector.names[i],
    ##      xlab="Lower tail index", ylab="Upper tail index");
    ## X <- seq(a, b, 0.01);
    ## lines(X, X, col="#FF0000");
    ## grid(nx=20);
    ## dev.off();

    pdf(sprintf("/tmp/%s_eigenvalue_ratios.pdf", sector.names[i]));
    company.num[i] <- length(I);
    # R.trfm <- rank.transform(R, I);
    E <- eigen(t(R[, I]) %*% R[, I], only.values=TRUE);
    l <- length(E$values);
    plot(1:(l-1), (E$values[2:l] / E$values[1:(l-1)]), type="b",
         xlab="", ylab="", xlim=c(1, 100), ylim=c(0.1, 1),
         main=sprintf("\"%s\"", sector.names[i]));
    grid(nx=20);
    dev.off();
}

pdf(sprintf("../papers/SP500_eigenvalue_ratios.pdf"));
E <- eigen(t(R) %*% R, only.values=TRUE);
l <- length(E$values);
plot(1:(l-1), (E$values[2:l] / E$values[1:(l-1)]), type="b",
     xlab="", ylab="", xlim=c(1, 100), ylim=c(0.1, 1));
dev.off();

################################################
## Remove the companies with small tail indices
################################################
min.index <- rep(NA, p);
for (i in 1:dim(tailIndices)[1]) {
    min.index[i] <- tailIndices[i,];
}

pdf("../papers/eigenvalue_ratios_selective.pdf");
colors <- rainbow(7);
largest <- rep(NA, 7);
E <- eigen(t(R) %*% R, only.values=TRUE);
largest[1] <- E$values[2]/E$values[1];
plot(1:(p-1), log10(E$values[2:p] / E$values[1:(p-1)]), type="b",
     xlab="", ylab="", xlim=c(1, 50), ylim=log10(c(0.1, 1)), col=colors[1]);
i <- 2;
thresholds <- seq(2, 2.5, by=0.1);
for (threshold in thresholds) {
    I <- which(min.index >= threshold);
    E <- eigen(t(R[, I]) %*% R[, I], only.values=TRUE);
    largest[i] <- E$values[2]/E$values[1];
    l <- length(E$values);

    ## pdf(sprintf("../papers/eigenvalue_ratios_%.1f.pdf", threshold));
    points(1:(l-1), log10(E$values[2:l] / E$values[1:(l-1)]), type="b", col=colors[i]);
    i <- i+1;
    ## dev.off();

    ## orders[i:(i+l-2)] <- 1:(l-1);
    ## ratios[i:(i+l-2)] <- E$values[2:l] / E$values[1:(l-1)];
    ## i <- i+l;
}
explanations <- c(
    "All included",
    expression(min*group("{", list(alpha["lower"], alpha["upper"]), "}") >= 2.0),
    expression(min*group("{", list(alpha["lower"], alpha["upper"]), "}") >= 2.1),
    expression(min*group("{", list(alpha["lower"], alpha["upper"]), "}") >= 2.2),
    expression(min*group("{", list(alpha["lower"], alpha["upper"]), "}") >= 2.3),
    expression(min*group("{", list(alpha["lower"], alpha["upper"]), "}") >= 2.4),
    expression(min*group("{", list(alpha["lower"], alpha["upper"]), "}") >= 2.5)
    );
legend("bottomright", legend=explanations, col=colors, lty=rep(1, 7));
dev.off();

## pdf(sprintf("../papers/eigenvalue_ratios_%.1f.pdf", threshold));
## dev.off();
## R.trfm <- rank.transform(R, 1:(dim(R)[2]));
## E <- eigen((n.records * p)^(-2) * t(R.trfm) %*% R.trfm, only.values=TRUE);
## lambda[1+length(sector.names),] <- E$values;


## grid(nx=20);
## legend("bottomright", legend=sector.names, col=colors, pch=seq(1, length(sector.names)));
## dev.off();

dbDisconnect(database);
