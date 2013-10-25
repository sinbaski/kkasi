clear all
close all
% Volvo B 2007-06-01 -- 2008-12-31 down
% Volvo B 2009-01-01 -- 2010-12-31 up
xi = 2.0512;

filename = 'data/ericsson_b.csv';
price = dlmread(filename, ',', [1, 4, 3560, 4]);
price = flipud(price); 


ret = price2ret(price);
r = (ret - mean(ret))/std(ret);

figure
[f, x] = ecdf(r);
plot(x, f, 'b');
R = randn(1e4, 1);
hold on
[f, x] = ecdf(R);
plot(x, f, 'r');
hold off
grid on
title('CDF of volvo daily ret. dist.');

figure
qqplot(r);
grid on

[p1, p2] = extract_xpnt(r);
fprintf('xpnt1=%f, xpnt2=%f\n', p1, p2);



% [a1, a2] = extract_xpnt(ret);
% fprintf('a1 = %e, a2 = %e\n', a1, a2);


% cen = linspace(min(ret), max(ret), 80);
% dr = cen(2) - cen(1);
% empirical = double(hist(ret, cen)) ./ (length(ret) * dr);
% log_empirical = log(empirical);

% pd = fitdist(ret, 'Normal');

% subplot(3, 1, 1);
% plot(cen, log_empirical, 'b.', cen, log(pdf(pd, cen)), 'r');
% title('return distribution log-log scale');
% grid on;

% subplot(3, 1, 2);
% plot(cen, empirical, 'b.', cen, pdf(pd, cen), 'r');
% title('return distribution');
% grid on;

% subplot(3, 1, 3);
% ret = ret ./ stdret;

% acf = autocorr(ret, 15);
% plot(0:15, acf, 'bo-');

