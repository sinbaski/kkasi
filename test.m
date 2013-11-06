clear all
close all

company = 'nordea_bank';
first_day = '2012-01-16';
last_day = '2012-04-20';
% dt in units of one minute
dt = 15;
% interval for calculating realized volatility. in seconds
delta = 30;

[ret, v] = get_intra_ret_simple(...
    company, first_day, last_day, dt, delta);
lv = log(v);
mu = mean(lv);
[d, nobs, tasy, sigasy, tols, sigols] = gph(lv - mu);
C = frct_cfct(d, 12);
y1 = frct_dfrc(lv - mu, C);
model = ar(y1, 3);
y2 = ar_filter(y1, model.A);
% [acf, x] = autocorr(ret.^2);
% autocorr(ret.^2);

% figure; qqplot(ret ./ sig);
% q = 5;
% p = 1;
% mdl0 = garch(p, q);
% [mdl, covariance, loglikelihood, info] = estimate(mdl0, ret(q+1:end-1), ...
%                                                   'E0', ret(1:q), ...
%                                                   'V0', h(q-p+1:q));

