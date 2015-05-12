clear all
close all
name='OMXS30';
T = 1000;
R = get_latest_data_from_index(name, T, true);
p = size(R, 2);

C = R' * R / T;
[portfolios, D] = eig(C);
X = R * portfolios(:, end-1:end);
Y = R(:, 8);
coef = X \ Y;
res = Y - X * coef;
autocorr(res, 40);

