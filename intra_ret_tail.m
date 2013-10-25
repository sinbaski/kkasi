clear all
close all
% Volvo B 2007-06-01 -- 2008-12-31 down
% Volvo B 2009-01-01 -- 2010-12-31 up
xi = 2.0512;

% ret = get_intra_ret_simple('ericsson_b', '2012-01-16', '2012-03-30', 50);
ret = get_intra_ret('ericsson_b', '2012-01-16', '2012-03-30', 45*60, ...
                                  1/60, 1, 0);
r = (ret - mean(ret))/std(ret);

% Now plot the pdf.
% n = length(r);
% x = linspace(-4, 4, n/10);
% y = hist(r, x);
% dx = x(2) - x(1);
% epdf = tsmovavg(y./(n*dx), 's', 3);
% epdf(2) = epdf(3)/2;
% epdf(1) = epdf(3)/4;
% plot(x, epdf);


% Now extract the exponent.
[p1, p2] = extract_xpnt(r);

% x1 = x(x < -xi);
% x2 = x(x > xi);
% y1 = epdf(x < -xi);
% y2 = epdf(x > xi);

% func = @(pa, arg) 1 ./ abs(arg).^pa(2);

% [p1, res1] = lsqcurvefit(func, [1, xi^2], x1, y1);
% [p2, res2] = lsqcurvefit(func, [1, xi^2], x2, y2);

% figure
% plot(x1, y1, 'b', x1, func(p1, x1), 'r');
% hold on
% plot(x2, y2, 'g', x2, func(p2, x2), 'm');
% grid on
% legend('left tail', sprintf('1/x^{%f}', p1(2)), 'right tail', ...
%        sprintf('1/x^{%f}', p2(2)), 'Location', 'South');
% title('nordea bank ret. dist. on tails.');

% hold off
