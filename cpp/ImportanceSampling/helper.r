rm(list = ls());

a1 <- 0.11
a2 <- 1.0e-8;
b1 <- 0.88;
c1 <- 0.011399856945298031;
c2 <- 0.82048971669444515;
c3 <- 0.66217902853870969;
delta <- 0.0053583136227694748;
p <- 0.99;
q <- 1;
rho <- 0.80795847750865057;
b <- 1.2573123595505618;

data <- read.table("nu_samples.txt", header=FALSE);

F <- list(x=c((a1 + rho)/(1 + a1), 1), eta=c(b1^2/rho, b));

nu.eta.pdf <- function(eta)
{
    pdf <- c1 / delta;
    pdf <- pdf * eta / q;
    pdf <- pdf / (eta * a1 - q * a2^2 + q * a2 * b1 + a2^2);
    pdf <- pdf * (F$x[2] - F$x[1]);
    return(pdf);
}

nu.eta.cdf <- function(eta)
{
    t1 <- a2 ^ 2
    t2 <- t1 * q
    t3 <- a2 * q
    t7 <- log(eta * a1 + b1 * t3 + t1 - t2)
    t8 <- q * t7
    t9 <- rho ^ 2
    t10 <- t1 * t9
    t13 <- b1 * a2 * t9
    t18 <- rho * t1
    t19 <- b1 ^ 2
    t20 <- a1 * t19
    t22 <- 0.1e1 / rho
    t24 <- log(-t22 * (-b1 * rho * t3 + rho * t2 - t18 - t20))
    t25 <- q * t24
    t34 <- t24 * t1
    t52 <- -b1 * q * rho * t24 * a2 + rho * a2 * b1 * t8 + a1 * t9 * eta - eta * rho * a1 + q * rho * t34 + t1 * t9 * t24 - t1 * t9 * t7 + t1 * rho * t7 - t19 * a1 * rho - t10 * t25 + t10 * t8 + t13 * t25 - t13 * t8 - t18 * t8 - rho * t34 + t20
    t60 <- a1 ^ 2
    t65 <- -t22 / t60 / q / (1 + a1) / delta * t52 * c1
    return (t65);
}

nu.eta.propdf <- function(eta)
{
    return(2 * eta / c3);
}

nu.eta.proqtl <- function(y)
{
}

nu.eta.procdf <- function(eta)
{
    return ((eta^2 - b1^4/rho^2) / c3);
}


## dx <- (F$eta[2] - F$eta[1])/1000;
## integral <- 0;
## for (x in seq(from=F$eta[1], to=F$eta[2], by=dx)) {
##     integral <- integral + nu.eta.propdf(x) * dx;
## }
F <- ecdf(data$V2);
X <- seq(from=b1^2 / rho, to=b, length.out=1001);
plot(X, F(X), type='l', );
lines(X, nu.eta.cdf(X), col="#00FF00");




