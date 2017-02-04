rm(list=ls());
source("libxxie.r");
## graphics.off();
## X11();
database = dbConnect(MySQL(), user='root', password='q1w2e3r4',
                     dbname='avanza', host="localhost");
Res <- dbSendQuery(database,
                   sprintf(
                       paste("select day,",
                             "closing from KLAC_US where day",
                             "between '2000-01-01' and '2015-01-01'"))
                   );
price <- fetch(Res, n=-1);
ratio <- price$closing[-length(price$closing)]/price$closing[-1];
I <- which(ratio < 0.51 | ratio > 1.9);


X <- diff(log(price$closing));


dbClearResult(Res);
dbDisconnect(database);

X11();
pdf("../papers/FX/MO_price.pdf");
plot(as.Date(price$day), price$closing, t="l",
     main="LLTC", xlab="time",
     ylab="price");
dev.off();


