clear all
close all

company = 'nordea_bank';
first_day = '2012-01-16';
last_day = '2012-04-20';
% dt in units of one minute
% dt = 15min, delta=30sec
% dt = 45min, delta=90sec,
% dt = 30min, delta=60sec
dt = 15;
% interval for calculating realized volatility. in seconds
delta = 30;

[ret, v] = get_intra_ret_simple(...
    company, first_day, last_day, dt, delta);
lv = log(v);
mu = mean(lv);
lv = lv - mu;

s = 33;
w = ts_difference(lv, s, 1);
w = ts_difference(w, 1, 1);
acf = autocorr(w, 120) .* var(w);
acf = acf(2:end);

% Nordea
theta = zeros(1, s+1);
theta(1) = 0.6263;
theta(s) = 0.5870;
theta(s+1) = theta(1)*theta(s);
aV = mean([acf(s-1), acf(s+1)])/theta(1)/theta(s);

% Assume Gaussian distribution to obtain
% a preliminary parameter estimation
% model = arima('D', 1, 'Seasonality', s, 'MALags', 1, 'SMALags', 1, ...
%               'Constant', 0);
% model.Distribution = struct('Name', 'Gaussian');
% model = estimate(model, lv);
% [y, yV] = infer(model, lv);
% theta = cell2mat(model.MA);

% Now infer the innovations
n = length(w);
y = NaN(n, 1);
done = 0;
for t = n:-1:1
    k = 1;
    x = 0;
    while 1
        m = floor(k/s);
        j = 0:m;
        coef = sum(theta(1).^(k - s.*j) .* theta(s).^j);
        if coef < 5.0e-2
            break;
        end
        if t - k <= 0
            done = 1;
            break;
        end
        x = x + w(t - k) * coef;
        k = k + 1;
    end
    if done
        break;
    end
    y(t) = w(t) + x;
end
y = y(t+1:n);



type = cmp_johnson_su(y);

% match the moments of the infered residuals to the corresponding
% formula of Johnson Su. so as to estimate the parameters a, b, c, m
func = @(x) -loglikelihood(x, w);

param = NaN(1, 6);
param([5,6]) = -theta([1, s]);
coef = fmincon(func, coef, [], [], [], [], -0.99, 0.99);


% [acf, x] = autocorr(ret.^2);
% autocorr(ret.^2);

%% FIT to a GARCH model
% q = 5;
% p = 1;
% mdl0 = garch(p, q);
% [mdl, covariance, loglikelihood, info] = estimate(mdl0, ret(q+1:end-1), ...
%                                                   'E0', ret(1:q), ...
%                                                   'V0', h(q-p+1:q));

