clear all
close all
L1 = importdata('DAX.txt', '\t');
L1 = L1(:, 2);
% subplot(2, 2, 1);
% qqplot(L1);
% grid on
% ylabel('Quantile of the DAX sample');

% subplot(2, 2, 2)
% sim_t = trnd(3, length(L1)*4, 1);
% qqplot(sim_t, L1);
% grid on
% xlabel('quantiles of t(3)');
% ylabel('Quantile of the DAX sample');

% subplot(2, 2, 3)
% sim_t = trnd(3.5, length(L1)*4, 1);
% qqplot(sim_t, L1);
% grid on
% xlabel('quantiles of t(3.5)');
% ylabel('Quantile of the DAX sample');

% subplot(2, 2, 4)
% sim_t = trnd(4, length(L1)*4, 1);
% qqplot(sim_t, L1);
% grid on
% xlabel('quantiles of t(4)');
% ylabel('Quantile of the DAX sample');

X = sort(-L1, 'ascend');
n = length(X);
Y = zeros(n, 1);
for k = 1 : n-1
    Y(k) = mean(X(k+1:end) - X(k));
end
plot(X, Y, '+');
xlabel('u');
ylabel('Mean Excess function e(u)');
I = X > 0.01;
m3 = skewness(X(I)-0.01);
xi = fsolve(@(x) 2*(1 + x)*sqrt(1 - 2*x) - (1 - 3*x)*m3, 1/4);
sig = sqrt(var(X(I) - 0.01)/(1 - xi^2)/(1 - 2*xi));
mu = mean(X(I) - 0.01) - sig/(1 - xi);

[param, pci] = mle(X(I)-0.01, 'distribution', 'gp');
[nlogL, acov] = gplike(param, X(I));


[y, x] = ecdf(X(I) - 0.01);
yi = gpcdf(x, param(1), param(2), 0);
plot(x, 1 - y, x, 1 - yi);
grid on
xlabel('x');
ylabel('1 - F(x+u)');
legend('empirical', 'fitted GPD');
