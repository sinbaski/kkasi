rm(list=ls());

r <- 0.05;
b <- 1;
p <- 0.5;
phi <- 0.5;

utility <- function(x) {
    return(log(x));
}

preference <- function(alpha, K, phi) {
    f <- function(x) {
        u <- log((1 - phi) * exp(r) + phi * exp(x));
        y <- u * (1 + b)/(K - x)^(alpha + 1);
        return(y);
    }
    y <- integrate(f, -Inf, 0);
    y <- alpha * K^alpha * p * y;
    return(y);
}

