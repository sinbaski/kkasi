rm(list=ls());
graphics.off();
## source("libxxie.r");
library("doParallel");

norm <- function(X) {
    return(max(X));
}

estimateLambda <- function(omega, a, b, theta)
{
    N <- 2000;
    K <- 2000;
    cl <- makeCluster(detectCores());
    registerDoParallel(cl);

    beta <- rep(NA, N);
    E <- matrix(runif(2*K), nrow=2, ncol=K);

    for (j in 1:N) {
        alpha <- rep(NA, K);
        Z2 <- rnorm(K)^2;
        A <- foreach(k=1:K, .combine='cbind') %dopar% {
            M <- matrix(NA, 2, 2);
            M[1, 1] <- a[1] * Z2[k] + b;
            M[1, 2] <- a[2];
            M[2, 1] <- Z2[k];
            M[2, 2] <- 0;
            M;
        }
        alpha <- foreach(k=1:K, .combine='c') %dopar% {
            norm(A[, (2*(k-1) + 1) : (2*k)] %*% E[, k])^theta
        }
        beta[j] <- sum(alpha);
        ## Draw r.v. k*
        Q <- cumsum(alpha/beta[j]);
        U <- runif(K);
        E <- foreach (k=1:K, .combine='cbind') %dopar% {    
            l <- min(which(Q > U[k]));
            V <- A[, (2*(k-1) + 1) : (2*k)] %*% E[, l];
            V / norm(V);
        }
    }
    stopCluster(cl);
    return(mean(log(beta/K)));
}

omega <- 1.0e-7;

## a <- c(0.11, 1.0e-8);
## b <- 0.88;
a <- c(0.6, 1.0e-3);
b <- 5.0e-3;

t1 <- Sys.time();
theta <- seq(from=1.1, to=2.5, by=0.1);
Lambda <- rep(NA, length(theta));
##for (i in 1:length(theta)) {
i <- 1;
Lambda[i] <- estimateLambda(omega, a, b, theta[i]);
##}
t2 <- Sys.time();
time.used <- difftime(t2, t1);

write.table(format(Lambda, digits=2), file="Lambda.dat", row.names=FALSE, col.names=FALSE, quote=FALSE);
print(time.used);


