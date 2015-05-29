clear all
close all
name='OMXS30_components';
T = 1000;
[R, syms] = get_latest_data_from_index(name, T, '2015-05-22', true, false);
p = size(R, 2);

acf = NaN(10*log10(T), p);
C = R' * R / T;
[portfolios, D] = eig(C);
portfolios = portfolios * diag(1./sum(abs(portfolios)));
factors = R * portfolios;
for k = 1 : p
    A = autocorr(factors(:, k), 10*log10(T));
    acf(:, k) = A(2:end);
end
ptfl = portfolios(:, 1);
ret = R * ptfl;
bic = NaN(5, 1);
for a = 1 : 5
    mdl0 = arima('ARLags', 1:a, 'MALags', 1:6,...
                 'Distribution', struct('Name', 't', 'DoF', NaN) ...
                 );
    [mdl, paramCov, loglik, info] = estimate(mdl0, ret);
    bic(a) = (a + 6 + 3) * log(T) - 2*loglik;
    
end
mdl0 = arima('ARLags', 1:3, 'MALags', 1:6,...
             'Distribution', struct('Name', 't', 'DoF', NaN) ...
             );
n = T/5;
[mdl, paramCov, loglik, info] = estimate(mdl0, ret(1:T-n));

predicted = NaN(n, 1);
for k = 1 : n
    predicted(k) = forecast(mdl, 1, 'Y0', ret(1:T-n+k-1));
end
realized = ret(end-n+1:end);
J = abs(predicted) > 0;
sum(predicted(J) .* realized(J) > 0) / sum(J)


% ptfl = portfolios(:, end);
% ptfl = ptfl ./ sum(abs(ptfl));
% I = find(abs(ptfl) > 0);
% s = sum(abs(ptfl(I)));
% ptfl = ptfl(I) ./ s;
% X = R(:, 1) - R(:, I) * ptfl;
% autocorr(X);

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

