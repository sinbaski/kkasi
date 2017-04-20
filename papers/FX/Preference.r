rm(list=ls());
library(plot3D)
library(parallel)

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

    y1 <- (1 + b) * integrate(function(x) {
        U <- utility(consumption(x, phi, r)) * dt(x, nu);
        return(U);
    }, -Inf, 0)$value;

    y2 <- (1 + b) * integrate(function(x) {
        U <- utility(consumption(x, phi, r)) * dt(x, nu);
        return(U);
    }, 0, q)$value;
    y2 <- y2 + integrate(function(x) {
        U <- utility(consumption(0, phi, r-x)) * dt(x, nu);
        return(U);
    }, q, Inf)$value;
    y2 <- y2 + integrate(function(x) {
        U <- (1/x^2 + 1/nu)^(-nu/2-1/2) * x^(-nu);
        return(U);
    }, q, Inf)$value;
    
    y3 <- b * pt(q, nu) * utility(consumption(q, phi, r));
    return(c(y1, y2, y3));
}

pareto.optimal.alloc <- function(alpha, alpha.r, K, K.r) {
    result <- optimize(
        function(phi, alpha, alpha.r, K, K.r) {
            y <- pareto.preference(phi, alpha, alpha.r, K, K.r);
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

scales <- seq(from=0.001, to=0.02, length.out=50);
indices <- seq(from=2, to=5, length.out=40);
phi.hat <- matrix(NA, nrow=length(indices),
                  ncol=length(scales));
U <- matrix(NA, nrow=length(indices),
            ncol=length(scales));
K.r <- 0.01;
alpha.r <- 2.0;
for (i in 1:length(indices)) {
    for (j in 1 : length(scales)) {
        phi <- pareto.optimal.alloc(indices[i], alpha.r,
                               scales[j], K.r);
        phi.hat[i, j] <- phi;
        y <- pareto.preference(phi, indices[i], alpha.r,
                               scales[j], K.r);
        U[i, j] <- y[1];
    }
}

## nu <- seq(from=1.5, to=5, length.out=200);
## phi.hat <- matrix(NA, nrow=length(nu), ncol=4);
## U <- matrix(NA, nrow=length(nu), ncol=4);
## for (k in 1:4) {
##     b <- 0.5 * (k - 1);
##     for (i in 1:length(nu)) {
##         phi.hat[i, k] <- t.optimal.alloc(nu[i]);
##         y <- t.preference(phi.hat[i, k], nu[i]);
##         U[i, k] <- y[1] + y[2] - y[3];
##     }
## }

pdf("phi_hat_b_t_power.pdf");
colors <- c("black", "green", "blue", "red");
plot(nu, phi.hat[, 1], type="l", lwd=2,
     ## xlim=c(1.5, 6),
     ylim=c(min(phi.hat), max(phi.hat)),
     xlab=expression(nu), col=colors[1],
     ## main=expression(list(u(C)==ln(C))),
     main=expression(list(u(C)==-frac(2, sqrt(C)))),
     ylab=expression(hat(phi)));
for (k in 2:4) {
    lines(nu, phi.hat[, k], col=colors[k], lwd=2);
}
dev.off();

pdf("U_b_t_power.pdf");
colors <- c("black", "green", "blue", "red");
plot(nu, U[, 1], type="l", lwd=2,
     ## xlim=c(1.5, 6),
     ylim=c(min(U), max(U)),
     xlab=expression(nu), col=colors[1],
     ## main=expression(list(u(C)==ln(C))),
     main=expression(list(u(C)==-frac(2, sqrt(C)))),
     ylab=expression(G(nu)));
for (k in 2:4) {
    lines(nu, U[, k], col=colors[k], lwd=2);
}
legend("topright", col=colors,
       lwd=rep(2, 4),
       legend=c(
           expression(b==0),
           expression(b==0.5),
           expression(b==1.0),
           expression(b==1.5)
       )
       );
dev.off();



pdf("phi_hat_pareto5e-1.pdf")
filled.contour(indices, scales, phi.hat, nlevels=60,
               xlab=expression(alpha), ylab=expression(K),
               color=terrain.colors,
               main=expression(u(x)==-frac(2, sqrt(x)))
               );
dev.off();

pdf("preference_pareto4.pdf");
filled.contour(indices, scales, U, nlevels=60,
               xlab=expression(alpha), ylab=expression(K),
               color=terrain.colors,
               main=expression(u(x)==-frac(1, 4*x^4))
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

