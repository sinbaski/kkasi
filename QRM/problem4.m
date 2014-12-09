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
I = X > 0.02;
