rm(list=ls());
library(plot3D)

consumption <- function(x, phi, r) {
    return((1 - phi) * exp(r) + phi * exp(x));
}

power.utility <- function(x) {
    return(-abs(x)^(-xi)/xi);
}

Pareto <- function(x, alpha, alpha.r, K, K.r) {
    if (x <= 0) {
        return(p * K^alpha/(K - x)^alpha);
    } else {
        return(1 - (1 - p)* (K.r / (K.r + x))^alpha.r);
    }
}


pareto.preference <- function(phi, alpha, alpha.r, K, K.r) {
    q <- Inf;
    if (phi > 0) {
        q <- log(exp(r) + (delta.v - exp(r))/phi);
    }

    y1 <- integrate(function(x) {
        U <- utility(consumption(x, phi, r)) *
            alpha * K^alpha * p /(K - x)^(alpha + 1);
        return(U * (1 + b));
    }, -Inf, 0)$value;

    y2 <- integrate(function(x) {
        U <- utility(consumption(x, phi, r)) *
            alpha.r * K.r^alpha.r * (1 - p) /(K.r + x)^(alpha.r + 1);
        U[x < q] <- U[x < q] * (1 + b);
        return(U);
    }, 0, Inf)$value;

    y3 <- b * Pareto(q, alpha, alpha.r, K, K.r) *
        utility(consumption(q, phi, r));

    return(c(y1, y2, y3));
}

t.preference <- function(phi, nu) {
    q <- Inf;
    if (phi > 0) {
        q <- log(exp(r) + (delta.v - exp(r))/phi);
    }

    y1 <- integrate(function(x) {
        U <- utility(consumption(x, phi, r)) * dt(x, nu);
        return(U * (1 + b));
    }, -Inf, 0)$value;

    y2 <- integrate(function(x) {
        U <- utility(consumption(x, phi, r)) * dt(x, nu);
        U[x < q] <- U[x < q] * (1 + b);
        return(U);
    }, 0, Inf)$value;

    y3 <- b * pt(q, nu) * utility(consumption(q, phi, r));
    return(c(y1, y2, y3));
}

pareto.optimal.alloc <- function(alpha, alpha.r, K, K.r) {
    result <- optimize(
        function(phi, alpha, alpha.r, K, K.r) {
            y <- preference(phi, alpha, alpha.r, K, K.r);
            return(y[1] + y[2] - y[3]);
        },
        interval=c(0, 1),
        alpha=alpha, alpha.r=alpha.r,
        K=K, K.r=K.r, maximum=TRUE
    );
    return(result$maximum);
}

t.optimal.alloc <- function(nu) {
    result <- optimize(
        function(phi, nu) {
            y <- t.preference(phi, nu);
            return(y[1] + y[2] - y[3]);
        },
        interval=c(0, 1), nu=nu,
        maximum=TRUE
    );
    return(result$maximum);
}

r <- 0.01;
# q <- r;
b <- 0.01;
## b <- 0;
p <- 0.5;
delta.v <- exp(r)*1.05;
xi <- 0.5;

## alpha.r <- 3.5;
## alpha <- 3;
## K.r <- 0.5;
## K <- 0.1;

## utility <- function(x) log(x);
utility <- power.utility;

## scales <- seq(from=0.01, to=1, length.out=50);
## indices <- seq(from=2, to=5, length.out=40);
## phi.hat <- matrix(NA, nrow=length(indices),
##                   ncol=length(scales));
## U <- matrix(NA, nrow=length(indices),
##             ncol=length(scales));
## for (i in 1:length(indices)) {
##     for (j in 1 : length(scales)) {
##         phi <- pareto.optimal.alloc(indices[i], indices[i],
##                                scales[j], scales[j]);
##         phi.hat[i, j] <- phi;
##         y <- pareto.preference(phi, indices[i], indices[i],
##                                scales[j], scales[j]);
##         U[i, j] <- y[1];
##     }
## }

nu <- seq(from=1.5, to=5, length.out=100);
phi.hat <- rep(NA, length(nu));
U <- rep(NA, length(nu));
for (i in 1:length(nu)) {
    phi.hat[i] <- t.optimal.alloc(nu[i]);
    y <- t.preference(phi.hat[i], nu[i]);
    U[i] <- y[1] + y[2] - y[3];
}

pdf("phi_hat_b_t5e-1.pdf");
plot(nu, phi.hat, type="p",pch=22,
     xlab=expression(nu),
     main=expression(list(u(C)==-frac(2, sqrt(C)), b == 0.01)),
     ylab=expression(hat(phi)));
grid(nx=20);
dev.off();


pdf("phi_hat_pareto5e-1.pdf")
filled.contour(indices, scales, phi.hat, nlevels=60,
               xlab=expression(alpha), ylab=expression(K),
               color=terrain.colors,
               main=expression(u(C)==-frac(2, sqrt(C)))
               );
dev.off();

pdf("preference_pareto5e-1.pdf");
filled.contour(indices, scales, U, nlevels=60,
               xlab=expression(alpha), ylab=expression(K),
               color=terrain.colors,
               main=expression(u(C)==-frac(2, sqrt(C)))
               );
dev.off();



## surf3D(M$x, M$y, phi.hat);

M <- matrix(NA, nrow=60, ncol=3);
portions <- seq(from=0.05, to=0.95, length.out=60)
for (i in 1:length(portions)) {
    M[i, ] <- preference(portions[i], alpha, alpha.r, K, K.r);
}
plot(portions, M[, 1] + M[, 2] - M[, 3],
     type="l",
     xlab=expression(phi), ylab="preference");
grid(nx=20);

result <- optimize(function(phi, alpha, alpha.r, K, K.r) {
    y <- preference(phi, alpha, alpha.r, K, K.r);
    return(y[1] + y[2] - y[3]);
    },
    interval=c(0, 1),
    alpha=alpha, alpha.r=alpha.r,
    K=K, K.r=K.r, maximum=TRUE);

