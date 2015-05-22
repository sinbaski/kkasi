clear all
close all
name='OMXS30';
T = 1000;
R = get_latest_data_from_index(name, T, true);
p = size(R, 2);

C = R' * R / T;
[portfolios, D] = eig(C);
X = R * portfolios(:, end-1);
% X = R * portfolios(:, 20);
% Y = ts_difference(X(:, 2), [10, 1; 1, 1]);
Y = ts_difference(X, [10, 1; 1, 1]);
% X = R * ones(p, 1);
autocorr(Y, 60);
% Y = R(:, 8);
% coef = X \ Y;
% res = Y - X * coef;
% autocorr(res, 40);
ylim([-0.8, 1]);

