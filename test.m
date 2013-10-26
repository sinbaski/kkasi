clear all
close all

company = 'nordea_bank';
first_day = '2012-01-16';
last_day = '2012-04-20';
% dt in units of one minute
dt = 15;

[ret, sig] = get_intra_ret_simple(company, first_day, last_day, dt);
h = sig.^2;

% q = 5;
% p = 1;
% mdl0 = garch(p, q);
% [mdl, covariance, loglikelihood, info] = estimate(mdl0, ret(q+1:end-1), ...
%                                                   'E0', ret(1:q), ...
%                                                   'V0', h(q-p+1:q));

