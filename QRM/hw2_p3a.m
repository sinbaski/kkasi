clear all
d = 50;
T = 50000;
kendal = 0.4;
rho = sin(pi * kendal / 2);
Sigma = rho * ones(d, d);
for m = 1 : d
    Sigma(m, m) = 1;
end
Z = mvnrnd(zeros(1, d), Sigma, T);
U = normcdf(Z);
X = tinv(U, 3)./(100 * sqrt(3));
L = (1 - exp(X)) * ones(d, 1) * 100;
VaR = quantile(L, 0.99);
ES = mean(L(L >= VaR));
% Y = trnd(3, T, 1);
% qqplot(X(:, 27), Y); xlim([-30, 30]); ylim([-30, 30]); grid on