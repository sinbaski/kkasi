qMyDist <- function(U, a) {
    X <- rep(NA, length(U));
    I1 <- which(U <= 1/4)
    U1 <- U[I1];
    I2 <- which(U > 1/4 & U <= 3/4);
    U2 <- U[I2];
    I3 <- which(U > 3/4);
    U3 <- U[I3];
    if (length(I1) > 0)
        X[I1] <- -1/(4^(1 + 1/a) * U1^(1/a));
    # X[I1] <- -(4^(a+1) * U1)^(-1/a);
    # X[I2] <- U2/a - (1 + 1/a)/4;
    if (length(I2) > 0)
        X[I2] <- U2 - 1/2;
    if (length(I3) > 0)
        X[I3] <- 1/(4^(1 + 1/a) * (1 - U3)^(1/a));
    # X[I3] <- 1/X[I3];
    return(X);
}

rMyDist <- function(n, a) {
    U <- runif(n);
    return (qMyDist(U, a));
}

q3Point <- function(U, p, loc) {
    X <- rep(NA, length(U));
    I1 <- U < p;
    I3 <- U > 1 - p;
    I2 <- rep(TRUE, length(U)) & !I1 & !I3;
    
    X[I1] <- -loc;
    X[I2] <- 0;
    X[I3] <- loc;
    return(X);
}

r3Point <- function(n, p, loc) {
    U <- runif(n);
    return (q3Point(U, p, loc));
}
    
rank.transform <- function(R, indices) {
    ranks <- matrix(NA, nrow=dim(R)[1], ncol=length(indices));
    for (j in 1:length(indices)) {
        ranks[, j] <- -1/log(rank(R[,indices[j]])/(dim(R)[1]+1));
    }
    return(ranks);
}
