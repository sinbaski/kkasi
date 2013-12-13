clear all
close all

company = 'volvo_b';
first_day = '2013-10-10';
last_day = '2013-12-04';
% dt in units of one minute
% dt = 15min, delta=30sec
% dt = 45min, delta=90sec,
% dt = 30min, delta=60sec
dt = 102;
% interval for calculating realized volatility. in seconds
delta = 204;

s = 4;
h = 0;

[ret_a, v_a] = get_intra_ret_simple(...
    company, first_day, last_day, dt, delta);
ret = ret_a(1:end-h);
lv_a = log(v_a);

v = v_a(1:end-h);
lv = lv_a(1:end-h);

clear model0;
model0 = arima('ARLags', [1:5], 'SARLags', [1:6]*4,...
               'Distribution', 'Gaussian');
load('volvo_2hr.mat', 'info');
info.options.MaxFunEvals = 1e4;
% [model0, VarCov, LogL, info] = ...
%     estimate(model0, lv);
[model0, VarCov, LogL, info] = ...
    estimate(model0, lv, 'options', info.options);
[Y, V] = infer(model0, lv);
save('volvo_2hr.mat', 'info');

save('/tmp/models.mat', 'model0');
