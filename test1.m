clear all
%%intraday returns
company = 'ericsson_b';
first_day = '2013-10-10';
last_day = '2014-01-29';
dt = 15;
ret = get_intra_ret(company, first_day, last_day, dt);
fprintf('%s skewness: %.4f\n', company, skewness(ret));

lag = 16;
x = linspace(min(ret)*0.9, max(ret), 130);
y = hist(ret, x) / length(ret) / (x(2) - x(1));
x = tsmovavg(x, 's', lag);
x = x(lag:end);
y = tsmovavg(y, 's', lag);
y = y(lag:end);
% y1 = pdf('Gamma', x, k, theta);
% plot(x, y, spec{(j-1)/4+1});
plot(x, y, 'b', x, normpdf(x, mean(ret), std(ret)), 'r');
