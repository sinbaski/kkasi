clear all
close all

company = 'volvo_b';
% first_day = '2012-01-16';
% last_day = '2012-04-20';
first_day = '2013-10-10';
last_day = '2014-04-04';

% dt in units of one minute
% dt = 15min, delta=30sec
% dt = 30min, delta=60sec
% dt = 45min, delta=90sec,
dt = 30;
% interval for calculating realized volatility. in seconds
delta = 120;

[ret_cmplt, v_cmplt] = get_intra_ret_simple(...
    company, first_day, last_day, dt, delta);
ret_cmplt = ret_cmplt - mean(ret_cmplt);
inno_cmplt = ret_cmplt ./ v_cmplt;
lv_cmplt = log(v_cmplt);
N = length(ret_cmplt);

ret = ret_cmplt(1:floor(N*0.8));
lv = lv_cmplt(1:floor(N*0.8));
lv1 = lv_cmplt(length(lv)+1:end);
inno = inno_cmplt(1:floor(N*0.8));
% lv1 = lv(1:s:end);

% mu = mean(lv);
% lv = lv - mu;
s = 16;

w = ts_difference(lv, s, 1);
w = ts_difference(w, 1, 1);
wm = min(w);
% [w, lam] = boxcox(w - wm + 1.0e-3);
acf = autocorr(w, 120);
acf = acf(2:end);

% Nordea


% theta = zeros(1, s+1);
% m = mean([acf(s-1), acf(s+1)]);
% theta(1) = fzero(@(x) (1+x^2)*m/acf(s) + x, 0);
% theta(s) = fzero(@(x) (1+x^2)*m/acf(1) + x, 0);
% theta(s+1) = -theta(1)*theta(s);
% y = ma_infer(w, theta(1), theta(s), s);


% Assume Gaussian distribution to obtain
% a preliminary parameter estimation
MALags = [1:4];
% model = arima('MALags', MALags, 'D', 1, 'Seasonality', s);
model = arima('MALags', MALags, 'SMALags', [1]*s, ...
              'D', 1, 'Seasonality', s);

% model.Distribution = struct('Name', 'T', 'DoF', NaN);
model.Distribution = struct('Name', 'Gaussian');
model = estimate(model, lv);
[y, V] = infer(model, lv);

lv2 = NaN(length(lv_cmplt) - length(lv), 1);
c = 1;
for k = length(lv):length(lv_cmplt) - 1
    lv2(c) = forecast(model, 1, 'Y0', lv_cmplt(1:k));
    c = c + 1;
end

%% Fit y_t to Johnson Su Distribution
% param = NaN(1, length(MALags));
% c = 1;
% for k = MALags
%     param(c) = model.MA{k};
%     c = c + 1;
% end
% mmt = [mean(y), var(y), skewness(y), kurtosis(y)];
% type = cmp_johnson_su(y);

% func = @(x) objective(x, model, [1, s, s+1], w);

% % theta & Theta must be in (-1, 1)
% lb = -ones(1, 6) .* Inf;
% lb(5:6) = [-0.99, -0.99];
% ub = ones(1, 6) .* Inf;
% ub(5:6) = [0.99, 0.99];

% % the mean of the innovations must be 0

% sltn = fmincon(func, param, [], [], [], [], lb, ub,...
%                @(x) nonlcon(x, model, MALags, w));
% c = 1;
% for k = MALags
%     model.MA{k} = sltn(c);
%     c = c + 1;
% end

% y1 = ma_infer(w, sltn(1), sltn(2), s);
% jsp = johnson_su_params([mean(y1), var(y1), skewness(y1), ...
%                     kurtosis(y1)]);






% ma = zeros(1, s+1);
% ma(1) = sltn(1);
% ma(s) = sltn(2);
% ma(s+1) = -prod(sltn(1:2));

% clear profile
% profile = nordea_15min_profile(lv, w, y1, ma, ...
%                                johnson_su_struct(jsp));

%% do h number of 1-step-ahead forecast
% if 0
% w_f = NaN(h, 1);
% lv_f = NaN(h, 1);
% m = length(lv);
% for k = 1:h
%     [w_f(k), lv_f(k)] = profile.forecast(1);
%     profile.update(lv_cmplt(m+k));
% end
% end

% if 0
% T = 1e3;
% N = 30;
% fold = 1e3;
% ev = NaN(1, 3e4);
% for k = 1:fold
%     M = profile.simulate_ret(T, N);
%     sigmas = diag(1./std(M));
%     % Normalize the simulated returns
%     M = M*sigmas;
%     C = M'*M/T;
%     ev((k-1)*N+1 : k*N) = eig(C);
% end
% hist(ev, 1e3);

% x = [-1:1e-3:1];
% rpdf = profile.ret_pdf(x, 1);

% [acf, x] = autocorr(ret.^2);
% autocorr(ret.^2);

%% FIT to a GARCH model
mdl0 = garch('ARCHLags', [1:3, 16], 'GARCHLags', [1, 16], 'Distribution', 'Gaussian');
[mdl, covariance, loglikelihood, info] = estimate(mdl0, ret);
c = 1;
v3 = NaN(length(lv2), 1);
U = v_cmplt.^2;
for k = length(lv):length(lv_cmplt) - 1
    v3(c) = forecast(mdl, 1, 'V0', U(1:k), 'Y0', ret_cmplt(1:k));
    c = c + 1;
end
lv3 = log(v3)./2;


%% Compare the differences of variance forecast
d2 = lv2 - lv1;
d3 = lv3 - lv1;
d4 = mean(lv) - lv1;
[u2, x2] = ecdf(lv2 - lv1);
[u3, x3] = ecdf(lv3- lv1);
[u4, x4] = ecdf(d4);

subplot(1, 2, 1);
plot((x2(x2 < 0)), (u2(x2<0)), 'b', ...
     (x3(x3 < 0)), (u2(x3<0)), 'g',...
     (x4(x4 < 0)), (u4(x4<0)), 'r',...
     'Linewidth', 2);
grid on
title('Cummulatibve Distribution of Under-estimates');
xlabel('x');
ylabel('P($\ln \sigma^F_t - \ln \hat{\sigma}_t < x$)', 'Interpreter', ...
       'Latex', 'Fontsize', 14);

subplot(1, 2, 2);
plot((x2(x2 > 0)), (1-u2(x2>0)), 'b', ...
     (x3(x3 > 0)), (1-u3(x3 > 0)), 'g', ...
     (x4(x4 > 0)), (1-u4(x4 > 0)), 'r', ...
     'Linewidth', 2);
grid on
title('Complementary Distribution of Over-estimates');
xlabel('x');
ylabel('P($\ln \sigma^F_t - \ln \hat{\sigma}_t > x$)', 'Interpreter', ...
       'Latex', 'Fontsize', 14);

%% Compare The quotient of variance forecast
% q2 = lv2 ./ lv1;
% q3 = lv2 ./ lv1;
% q4 = mean(lv) ./ lv1;

% [u2, x2] = ecdf(q2);
% [u3, x3] = ecdf(q3);
% [u4, x4] = ecdf(q4);
% plot(log10(x2), log10(u2), log10(x3), log10(u3), log10(x4), log10(u4));
