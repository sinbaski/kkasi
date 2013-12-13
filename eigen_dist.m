clear all

T = 1e3;
N = 50;
fold = 4e3;

q = N/T;
ev = NaN(1, N*fold);
% mp_max = (1 + sqrt(q))^2;

%% eigen values generated from the volvo stochastic log-volatility
%% model.

spec = ['b', 'c', 'g', 'm', 'r', 'k'];
for i = [1:5]
    load(sprintf('./eigen-%d.mat', i));
    ev = reshape(ev, N, fold);
    maxeig = max(ev);
    [y, x] = ecdf(maxeig);
    y1 = diff(y) ./ diff(x);
    x1 = (x(1:end-1) + x(2:end))/2;
    plot(log(x1), log(y1), spec(i));
    hold on
end
hold off
grid on


%% White Wishart matrix eigen values
% for i = 1:fold
%     A = randn(N, T);
%     C = A*A';
%     ev(N*(i-1)+1 : N*i) = eig(C);
% end

% ev = reshape(ev, N, fold);
% maxeig = max(ev);
% p = RealWishartCDF(x, N, T);
% plot(x, y, x, p);

%% Approximation of the Tracy-Widom with a Gamma dist.
% a1 = -1/2;
% a2 = a1;
% n = N;
% p = T;
% mu = (sqrt(n+a1) + sqrt(p+a2))^2;
% sigma = sqrt(mu) * (1/sqrt(n+a1) + 1/sqrt(p+a2))^(1/3);
% maxeig1 = (maxeig - mu)/sigma;

% moments = [mean(maxeig1), var(maxeig1), skewness(maxeig1)];
% k = 4/moments(3)^2;
% theta = prod(moments(2:3))/2;
% alpha = k*theta - moments(1);

% maxeig1 = maxeig1 + alpha;
% % x = linspace(min(maxeig1), max(maxeig1), 60);
% % y = hist(maxeig1, x) / length(maxeig1) / (x(2) - x(1));

% [y, x] = ecdf(maxeig1);

% y1 = cdf('Gamma', x, k, theta);
% % plot(x, y);

% %% Fit data to Gamma dist.
% m = min(maxeig1);
% pd = fitdist(maxeig1' - m, 'Gamma');
% y2 = cdf(pd, x-m);

% plot(log(x-m), log(y), log(x-m), log(y2));
