rm(list=ls());

consumption <- function(x, phi, r) {
    return((1 - phi) * exp(r) + phi * exp(x));
}

utility <- function(x) {
    return(log(x));
    ## abs(x)^0.5/(1 - 0.5);
    ## -abs(x)^(-3)/3;
}

F <- function(x, alpha, alpha.r, K, K.r) {
    if (x <= 0) {
        return(p * K^alpha/(K - x)^alpha);
    } else {
        return(1 - (1 - p)* (K.r / (K.r + x))^alpha.r);
    }
}

## density.X <- function(x, r, phi, b) {
##     x[x < 0]
## }

## q <- function(r, phi) {
##     return(log(exp(r) + (delta.v - exp(r))/phi));
## }

## f <- function(x, K, alpha, b) {
##     y <- alpha * K^alpha * /(K - x)^(alpha + 1);
##     return(y);
## }

## g <- function(x, K.r, alpha.r) {
##     u <- utility((1 - phi) * exp(r) + phi * exp(x));
##     I <- which(x <= q);
##     u[I] <- u[I] * (1 + b);
##     u <- u / (K.r + x)^(alpha.r + 1);
##     return(u);
## }

preference <- function(phi, alpha, alpha.r, K, K.r) {
    if (phi > 0) {
        q <- log(exp(r) + (delta.v - exp(r))/phi);
    } else {
        q <- Inf;
    }

    y1 <- integrate(function(x) {
        U <- utility(consumption(x, phi, r)) * alpha * K^alpha * p /(K - x)^(alpha + 1);
        return(U * (1 + b));
    }, -Inf, 0)$value;

    y2 <- integrate(function(x) {
        U <- utility(consumption(x, phi, r)) * alpha.r * K.r^alpha.r * (1 - p) /(K.r + x)^(alpha.r + 1);
        U[x < q] <- U[x < q] * (1 + b);
        return(U);
    }, 0, Inf)$value;

    ## return(c(y1, y2, y1 + y2 - F(q, alpha, alpha.r, K, K.r) * b * utility(delta.v)));
    return(y1 + y2 - F(q, alpha, alpha.r, K, K.r) * b * utility(delta.v));
}

r <- 0.02;
b <- 0.01;
p <- 0.5;
phi <- 0.8;
delta.v <- exp(r)*1.05;
alpha.r <- 3.5;
alpha <- 3.5;
K.r <- 0.001;
K <- 0.001;

M <- matrix(NA, nrow=20, ncol=3);
counter <- 1;
for (phi in seq(from=0, to=1, length.out=20)) {
    M[counter, ] <- preference(phi, 3, 1.8, 0.03, 0.03);
    counter <- counter + 1;
}
plot(M[, 2]);
result <- optimize(preference, interval=c(0, 1), alpha=2, alpha.r=3.5, K=2.0e-6, K.r=2.0e-6, maximum=TRUE);
## result <- optimize(function(x) -(x-1)^2 + 1, c(0, 2), maximum=TRUE);
result
