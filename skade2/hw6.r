rm(list = ls());

data = read.csv("bulkclaim.csv", header=TRUE, sep=",", dec=".");
data = as.matrix(data);
loss = log(data[, 2]);
Ded = log(data[, 3]);

n = length(loss);
BIC = rep(NA, 7);

c = 1;
for (m in (1:3)) {
   S = combn((4:6), m);
   for (j in (1:dim(S)[2])) {
        Z = data[, S[, j]];
        p = dim(S)[1];
        fun1 <- function(theta) {
            l = length(theta);
            intercept = theta[1];
            tau = theta[l];
            bita = theta[2:(l-1)];

            if (p >= 2) {
                A = loss - Z %*% bita - intercept;
            } else {
                A = loss - Z * bita - intercept;
            }
            loglikelihood = -n/2*log(tau) - 1/(2*tau) * sum(A^2);
            ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.BaseNamespaceEnv[['.ESSBP.']][["@18@"]]))##:ess-bp-end:##
            accum = 0;
            if (p >= 2) {
                for (k in (1:n)) {
                    accum = accum + log(pnorm((Z[k, ] %*% bita - Ded[k])/sqrt(tau)));
                }
            } else {
                for (k in (1:n)) {
                    accum = accum + log(pnorm(sum((Z[k])) * bita - Ded[k])/sqrt(tau));
                }
            }
            loglikelihood = loglikelihood - accum;
            -loglikelihood;
        }


        bita0 = solve(t(Z)%*%Z) %*% t(Z) %*% loss;
        tau0 = sum((loss - Z %*% bita0)^2) / n;
        v = -fun1(c(0, bita0, tau0));
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.BaseNamespaceEnv[['.ESSBP.']][["@17@"]]))##:ess-bp-end:##
        if (v < Inf) {
            ret1 = optim(c(0, bita0, tau0), fun1, hessian=TRUE);

            par = ret1$par;
            sigma2Hat = sum((loss - Z %*% par[2:2+p-1] - par[1])^2)/n;
            BIC[c] = n*log(sigma2Hat) + p*log(n);
        } else {
            BIC[c] = Inf;
        }
        c = c + 1;
    }
}
