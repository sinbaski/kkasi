clear all
close all
T = 120;
horizon = 3;
day = '2014-07-18';
mysql = get_mysql();
stmt = sprintf(['select closing from DAX order by day;']);
data = fetch(mysql, stmt);
close(mysql);

% T = size(data, 1);
S = cell2mat(data(:, 1));
% price = flipud(price);
r = price2ret(S);
avg_span = 5;
forecasts = NaN(length(S) - (horizon - 1) - T, horizon);
for d = 1: size(forecasts, 1)
    price = S(d : d + T -1);
    %% Analysis based on auto-correlations in the returns
    % s = 14;
    % w = ts_difference(r, s, 1);

    % figure;
    % subplot(2, 1, 1);
    % autocorr(r, floor(length(r)/2));
    % subplot(2, 1, 2);
    % autocorr(w, floor(length(w)/2));

    % model = arima('Seasonality', s, 'SMALags', [1]*s, ...
    %              'Distribution', 'Gaussian');
    % model = estimate(model, r);
    % [res, V] = infer(model, r);
    % [y1, ymse1] = forecast(model, horizon);

    %% Analysis based on correlation between the return and
    % the price level
    %avg_span = [5, 20, 60];
    %avg_span = 10:5:60;
    avg = NaN(T, length(avg_span));
    B = NaN(length(avg_span), 1);
    for n = 1 : length(avg_span)
        l = T - avg_span(n) + 1;
        avg(:, n) = log(price) - log(tsmovavg(price', 's', avg_span(n))');
    end
    X = avg(max(avg_span):end, :);
    A = cov(X);
    [V, D] = eig(A);
    Factors = X * V;
    for n = 1 : length(avg_span)
        m = size(X, 1);
        U = cov(r(end-m+1:end), Factors(:, n));
        B(n) = U(1, 2);
    end
    Theta = inv(D) * B;
    y2 = NaN(horizon, 1);
    y = NaN(horizon, 1);
    u = NaN(1, length(avg_span));
    for n = 1: horizon
        y2(n) = Factors(end, :) * Theta;
        y(n) = y2(n);
        price(end + 1) = price(end) * (1 + y(n));
        for m = 1: length(avg_span)
            u(m) = log(price(end)) - ...
                   log(mean(price(end-avg_span(m)+1 : end)));
        end
        X(end+1, :) = u;
        Factors(end+1, :) = u * V;
    end
    forecasts(d, :) = y';
end
difference = NaN(size(forecasts, 1), horizon);
for n = 1 : horizon
    difference(:, n) = forecasts(:, n) - r(T+n-1:end-horizon+n);
end
epdf(difference(:, 1), 1, min(difference(:, 1)),...
     max(difference(:, 1)), 100, 'g-');
grid on
