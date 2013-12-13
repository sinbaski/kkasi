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
s = 33;
h = 0;

[ret_cmplt, v_cmplt] = get_intra_ret_simple(...
    company, first_day, last_day, dt, delta);
ret = ret_cmplt(1:end-h);
lv_cmplt = log(v_cmplt);

v = v_cmplt(1:end-h);
lv = lv_cmplt(1:end-h);

lv1 = lv(1:s:end);

% mu = mean(lv);
% lv = lv - mu;


w = ts_difference(lv, s, 1);
w = ts_difference(w, 1, 1);
acf = autocorr(w, 120);
acf = acf(2:end);

% Nordea

if 1
    theta = zeros(1, s+1);
    m = mean([acf(s-1), acf(s+1)]);
    theta(1) = fzero(@(x) (1+x^2)*m/acf(s) + x, 0);
    theta(s) = fzero(@(x) (1+x^2)*m/acf(1) + x, 0);
    theta(s+1) = -theta(1)*theta(s);
    y = ma_infer(w, theta(1), theta(s), s);
else

    % Assume Gaussian distribution to obtain
    % a preliminary parameter estimation
    model = arima('MALags', [1], 'SMALags', [1]*33, ...
                  'ARLags', [1,2], 'SARLags', [1,2]*33);
    % model = arima('MALags', [1], 'SMALags', [1]*33, ...
    %               'D', 1, 'Seasonality', 33);

    model.Distribution = struct('Name', 'T', 'DoF', NaN);
    model = estimate(model, lv);
    [Y, V] = infer(model, lv);
    % theta = cell2mat(model.MA);
end

param = NaN(1, 2);
param(1) = theta(1);
param(2) = theta(s);

mmt = [mean(y), var(y), skewness(y), kurtosis(y)];
type = cmp_johnson_su(y);

func = @(x) objective(x, w);

% theta & Theta must be in (-1, 1)
lb = -ones(1, 6) .* Inf;
lb(5:6) = [-0.99, -0.99];
ub = ones(1, 6) .* Inf;
ub(5:6) = [0.99, 0.99];

% the mean of the innovations must be 0

sltn = fmincon(func, param, [], [], [], [], lb, ub,...
               @(x) nonlcon(x, w));
y1 = ma_infer(w, sltn(1), sltn(2), s);
jsp = johnson_su_params([mean(y1), var(y1), skewness(y1), ...
                    kurtosis(y1)]);
ma = zeros(1, s+1);
ma(1) = sltn(1);
ma(s) = sltn(2);
ma(s+1) = -prod(sltn(1:2));

clear profile
profile = nordea_15min_profile(lv, w, y1, ma, ...
                               johnson_su_struct(jsp));

%% do h number of 1-step-ahead forecast
if 0
w_f = NaN(h, 1);
lv_f = NaN(h, 1);
m = length(lv);
for k = 1:h
    [w_f(k), lv_f(k)] = profile.forecast(1);
    profile.update(lv_cmplt(m+k));
end
end

if 0
T = 1e3;
N = 30;
fold = 1e3;
ev = NaN(1, 3e4);
for k = 1:fold
    M = profile.simulate_ret(T, N);
    sigmas = diag(1./std(M));
    % Normalize the simulated returns
    M = M*sigmas;
    C = M'*M/T;
    ev((k-1)*N+1 : k*N) = eig(C);
end
hist(ev, 1e3);

hold on;

x = linspace(min(ev), max(ev), 1e3);
rho = MarcenkoPasturPDF([N/T, 1], x);
plot(x, rho*(x(2) - x(1))*length(ev), 'r');
end



% x = [-1:1e-3:1];
% rpdf = profile.ret_pdf(x, 1);

% [acf, x] = autocorr(ret.^2);
% autocorr(ret.^2);

%% FIT to a GARCH model
% q = 5;
% p = 1;
% mdl0 = garch(p, q);
% [mdl, covariance, loglikelihood, info] = estimate(mdl0, ret(q+1:end-1), ...
%                                                   'E0', ret(1:q), ...
%                                                   'V0', h(q-p+1:q));

