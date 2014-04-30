function [p1, p2] = extract_xpnt(r, lb, ub)
%% r MUST BE NORMALIZED!!!!

% Volvo B 2007-06-01 -- 2008-12-31 down
% Volvo B 2009-01-01 -- 2010-12-31 up
% xi = 2.0512;

% Now extract the exponent.

[y, x] = ecdf(r);
x = x(2:end-1);
y = y(2:end-1);
x1 = x(x < -lb & x > -ub);
x2 = x(x > lb & x < ub);
y1 = y(x < -lb & x > -ub);
y2 = y(x > lb & x < ub);

% func = @(pa, arg) (1 - pa(1))*arg - pa(2);
% [P1, res1] = lsqcurvefit(func, [lb^2, log(1 - lb^2)], log(-x1), log(y1));
% [P1, res1] = lsqcurvefit(func, [lb^2, log(1 - lb^2)], log(x2), log(y2));

P1 = polyfit(log(-x1), log(y1), 1);
P2 = polyfit(log(x2), log(1-y2), 1);
p1 = -P1(1);
p2 = -P2(1);

figure
plot(log(-x1), log(y1), 'b', log(x2), log(1-y2), 'g');
hold on
plot(log(-x1), polyval(P1, log(-x1)), 'r');
plot(log(x2), polyval(P2, log(x2)), 'm');
hold off
grid on
legend(...
sprintf('ln(P(r < x)) vs. ln(-x)  x<-%.2e', lb),...
sprintf('ln(P(r > x)) vs. ln(x)  x>%.2e', lb),...
sprintf('Linear Fit'),...
sprintf('Linear Fit'),...
'Location', 'Southwest');
title(sprintf('left tail: 1/x^{%3.2f} --- right tail: 1/x^{%3.2f}', ...
              p1, p2));
xlabel('x');
