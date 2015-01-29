clear all
close all
q = 1;
sig = 0.5;
phi = 0.5;
v = sig^2/(1 - phi^2);
% variance(m, n) = v;
load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
              'sig%.4f-phi%.4f.mat'], q, sig, phi), 'ev');

% ev = reshape(ev, 1, prod(size(ev)));
eigmax = sort(max(ev), 'ascend')';
n = length(eigmax);
meanExcess = NaN(length(eigmax) - 1, 2);
H = NaN(length(eigmax) - 1, 1);
for k = 1 : length(eigmax)-1
    meanExcess(k, 2) = mean(eigmax(k+1:end) - eigmax(k));
    meanExcess(k, 1) = eigmax(k);
    H(n - k) = 1/mean(log(eigmax(k+1:end) ./ eigmax(k)));
end

subplot(1, 2, 1);
plot(meanExcess(:, 1), meanExcess(:, 2), 'x');
grid on
xlabel('u: threshold');
ylabel('Mean Excess Function e_F(u)');

subplot(1, 2, 2);
plot(1 : length(eigmax)-1, H);
xlim([1, 400]);
grid on
xlabel('#upper-order statistics');
ylabel('Hill Estimator');

gp = mle(eigmax(eigmax > 38), 'distribution', 'gp');

I = meanExcess(:, 1) > 38;
P = polyfit(meanExcess(I, 1), meanExcess(I, 2), 1);
Q = P ./ (1 + P(1));

params = mle(eigmax(eigmax > 38), 'pdf', @tmp, ...
             'start', Q, 'lowerbound', [0, 0]);

