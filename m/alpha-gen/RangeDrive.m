clear all
close all
mysql = get_mysql();
data = cell2mat(fetch(mysql, ['select low, closing, high from DAX order by ' ...
                    'day desc limit 250;']));
data = flipud(data);
close(mysql);
ranges = log(data(:, 3)) - log(data(:, 1));
closing = log(data(:, 2));
ret = price2ret(data(:, 2));

%% Fit a ARMA model to the ranges
model_init = arima('ARLags', 1:3, 'Distribution', 'Gaussian');
model = estimate(model_init, ranges);
[residuals, condvar] = infer(model, ranges);
Z = residuals ./ sqrt(condvar);

numlags = 10;
coef = NaN(numlags, 1);
lags = NaN(numlags, 1);
counter = 0;
for k = 1 : numlags
    % c = cov(Z(1:end-k), closing(k+1 : end) - closing(k : end-1));
    if abs(corr(Z(1:end-k), ret(k:end))) > 5.0e-2
        M = cov(Z(1:end-k), ret(k:end));
        coef(k) = M(1, 2);
        counter = counter + 1;
        lags(counter) = k;
    else
        break;
    end
end
lags = lags(~isnan(lags));
coef = coef(~isnan(coef));

deviations = NaN(129, 1);
for k = 121 : length(ret)
    deviations(k-120) = coef' * Z(k-lags+1) - ret(k);
end
