clear all
close all
name='OMXS30_components';
T = 1000;
[R, syms] = get_latest_data_from_index(name, T, '2015-05-22', true, false);
p = size(R, 2);

acf = NaN(60, p);
C = R' * R / T;
[portfolios, D] = eig(C);
portfolios = portfolios * diag(1./sum(abs(portfolios)));
factors = R * portfolios;
for k = 1 : p
    A = autocorr(factors(:, k), 60);
    acf(:, k) = A(2:end);
end

ptfl = portfolios(:, end);
ptfl = ptfl ./ sum(abs(ptfl));
I = find(abs(ptfl) > 0);
s = sum(abs(ptfl(I)));
ptfl = ptfl(I) ./ s;
X = R(:, 1) - R(:, I) * ptfl;
autocorr(X);

% Y = X;
% X = R * portfolios(:, 20);
% Y = ts_difference(X(:, 2), [10, 1; 1, 1]);
% Y = ts_difference(X, [10, 1; 1, 1]);
% X = R * ones(p, 1);
% autocorr(Y, 60);
% acf = autocorr(Y, 60);
% Y = R(:, 8);
% coef = X \ Y;
% res = Y - X * coef;
% autocorr(res, 40);
% ylim([-0.8, 1]);

