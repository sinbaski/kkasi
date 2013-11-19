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

if 1
    theta = zeros(1, s+1);
    % theta(1) = 0.6263;
    % theta(s) = 0.5870;
    theta(1) = 0.845552;
    theta(s) = 0.804139;
    % theta(1) = 0.878236;
    % theta(s) = 0.912821;
    theta(s+1) = theta(1)*theta(s);
    theta(s-1) = theta(s+1);
    V = mean([acf(s-1), acf(s+1)])/theta(1)/theta(s);
    
    y = ma_infer(w, theta(1), theta(s), s);

else

    % Assume Gaussian distribution to obtain
    % a preliminary parameter estimation
    model = arima('MALags', [1], 'SMALags', [33], ...
                  'Constant', 0);
    model.Distribution = struct('Name', 'Gaussian');
    model = estimate(model, w);
    [y, yV] = infer(model, w);
    % theta = cell2mat(model.MA);
end

param = NaN(1, 6);
param(5) = model.MA{1};
param(6) = model.SMA{33};
mmt = [mean(y), var(y), skewness(y), kurtosis(y)];

func = @(x) johnson_su_moments34(x) - mmt(3:4);
param([1:2]) = fsolve(func, [-0.1922, 1.445]);
m = johnson_su_moments12([param([1:2]), 1, 0]);
param(3) = sqrt(mmt(2)/m(2));
param(4) = mmt(1) - m(1)*param(3);
type = cmp_johnson_su(y);

func = @(x) objective(x, w);

lb = -ones(1, 6) .* Inf;
lb(5:6) = [-0.99, -0.99];
ub = ones(1, 6) .* Inf;
ub(5:6) = [0.99, 0.99];

% options = optimoptions('fmincon','GradObj','on');
% problem = struct('objective', func, 'x0', param, 'lb', lb, 'ub', ub, ...
%                  'options', options);
sltn = fmincon(func, param, [], [], [], [], lb, ub, []);


% [acf, x] = autocorr(ret.^2);
% autocorr(ret.^2);

%% FIT to a GARCH model
% q = 5;
% p = 1;
% mdl0 = garch(p, q);
% [mdl, covariance, loglikelihood, info] = estimate(mdl0, ret(q+1:end-1), ...
%                                                   'E0', ret(1:q), ...
%                                                   'V0', h(q-p+1:q));

