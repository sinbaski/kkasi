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

[r, v] = get_intra_ret_simple(company, first_day, last_day, dt, delta);
ep = r - mean(r);
z = ep ./ v;
probplot(log(v));
