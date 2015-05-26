clear all
close all
name='OMXS30';
T = 1000;
R = get_latest_data_from_index(name, T, true);
p = size(R, 2);

acf = NaN(60, p);
C = R' * R / T;
[portfolios, D] = eig(C);
for k = 0 : p-1
    X = R * portfolios(:, end-k);
    A = autocorr(X, 60);
    acf(:,k+1) = A(2:end);
end
I = abs(portfolios(:, 1)) > 0.05;
ptfl = zeros(p, 1);
ptfl(I) = portfolios(I, 1);
X = R * ptfl;
% X1 = ts_difference(X, [5, 2; 1, 1]);
% autocorr(X1, 40);
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

