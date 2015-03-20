rm(list = ls());

data = read.csv("bulkclaim.csv", header=TRUE, sep=",", dec=".");
claim = (data[, 2] - data[, 3]) * data[, 4] * 1.0e-6;

## loglikelihood of lognormal
X = log(claim);
fun1 <- function(theta) {
    A = dlnorm(claim, meanlog=theta[1], sdlog=theta[2]);
    -sum(log(A));
}
ret1 = optim(c(mean(X), sd(X)), fun1, hessian=TRUE);
kstest <- ks.test((log(claim)-ret1$par[1])/ret1$par[2], "pnorm");

a = mean(claim);
b = var(claim);

## loglikelihood of gamma
## if (a^2 / b < 1) {
##     sprintf("Cannot be Gamma distributed. a = %.2f\n", a^2/b);
##     ret2 <- list(value=Inf);
## } else {
##     fun2 <- function(theta) {
##         A = dgamma(claim, shape=theta[1], rate=theta[2]);
##         -sum(log(A));
##     }
##     ret2 <- optim(c(max(1, a^2/b), a/b), fun2, hessian=TRUE);
## }


## Weibull distribution
f <- function(x) (gamma(1 + 2/x)/gamma(1+1/x)^2 - 1 - b/a^2);
shape = uniroot(f, c(0.1, 10))$root;
scale = a/gamma(1 + 1/shape);
fun3 <- function(theta) {
    A = dweibull(claim, shape=theta[1], scale=theta[2]);
    -sum(log(A));
}
ret3 <- optim(c(shape, scale), fun3, hessian=TRUE);
N1 = rpois(

X = pweibull(claim, shape=ret3$par[1], scale=ret3$par[2]);
kstest <- ks.test(qnorm(X), "pnorm");

## Fit to Pareto
X = sort(claim);
l = length(X);
ef = rep(NA, l-1);
for (n in (1:l-1)) {
    ef[n] <- mean(X[(n+1):l] - X[n]);
}
plot(X[2:l], ef, type='p');
dpareto <- function(x, exponent=1, scale=1) {
    
}
