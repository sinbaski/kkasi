rm(list=ls());
graphics.off();
## require("tseries");
require("fGarch");
require("mvtnorm");
source("libxxie.r");

currencies <- c(
    "AUD_SEK_Rates",
    "CAD_SEK_Rates",
    "CHF_SEK_Rates",

    "CNY_SEK_Rates",
    "CZK_SEK_Rates",
    "DKK_SEK_Rates",

    "EUR_SEK_Rates",
    "GBP_SEK_Rates",
    "HKD_SEK_Rates",

    "HUF_SEK_Rates",
    "JPY_SEK_Rates",
    "KRW_SEK_Rates",

    "MAD_SEK_Rates",
    "MXN_SEK_Rates",
    "NOK_SEK_Rates",

    "NZD_SEK_Rates",
    "SGD_SEK_Rates",
    "USD_SEK_Rates"
    );


database = dbConnect(MySQL(), user='root', password='q1w2e3r4',
                     dbname='avanza', host="localhost");
i <- 4;
results <- dbSendQuery(database,
                       sprintf("select rate from %s order by day;",
                               currencies[i]));
X <- fetch(results, n=-1)[[1]];
dbClearResult(results);
X <- diff(X);
H <- hillPlot(X^2);

H$alpha[100];

dbDisconnect(database);
