function code = log_log_plot(returns)

N = 100;
returns = (returns - mean(returns)) ./ std(returns);
cen = linspace(min(returns), max(returns), N);

Y = hist(returns, cen);
% indices = (Y > 1);
% Y = double(Y(indices));
% X = cen(indices);
dr = cen(2) - cen(1);

use_movavg = 0;
if use_movavg
    Y = tsmovavg(Y, 's', 5, 2);
    Y = Y(~isnan(Y));
    X = tsmovavg(X, 's', 5, 2);
    X = X(~isnan(X));
    dr = X(2) - X(1);
end

tot = sum(Y);

empirical = Y ./ (tot * dr);
left_tail = (X < 0);
right_tail = (X > 0);

X2 = linspace(0, max(abs(X)), 60);
Y2 = pdf('Normal', X2, 0, 1);

x0 = 2.0512;

func = @(p, x) 1./abs(x).^norm(p);
[alpha1, resnorm] = lsqcurvefit(func, x0^2, ...
                                -X(X < -x0), empirical(X < -x0));
alpha1 = norm(alpha1);

[alpha2, resnorm] = lsqcurvefit(func, x0^2, ...
                                X(X > x0), empirical(X > x0));
alpha2 = norm(alpha2);

xl = min(X);
xr = max(X);
L = min(X):dr:-x0;
R = x0:dr:max(X);

hdl = figure;
plot(log(-X(left_tail)), log(empirical(left_tail)), 'bx', ...
     log(X(right_tail)), log(empirical(right_tail)), 'g+', ...
     log(X2), log(Y2), 'r-', ...
     log(abs(L)), -alpha1*log(abs(L)), 'b', ...
     log(R), -alpha2*log(R), 'g');
grid on;
stmt = sprintf(['%s %dmin Return distribution %s -- %s log-log scale'], ...
               strrep(company, '_', ' '), interval, start_day, ...
               end_day);
title(stmt);
legend('left tail', 'right tail', ...
       'standard Gaussian', ...
       sprintf('power-law with xpnt %.4f', alpha1), ...
       sprintf('power-law with xpnt %.4f', alpha2), ...
       'Location', 'SouthWest');
saveas(hdl, sprintf('pics/%s_%dmin_ret_%s-%s_log.pdf', company, ...
                    interval, start_day, end_day));

% To separate the distribution into a central Gaussian and two
% power-law tails 1/x^a, the separation point is at
% x^2(1 - ln x^2) + ln(2pi) = 0
hdl = figure;
plot(X, empirical, 'bx', ...
     L, 1./abs(L).^alpha1, 'b', ...
     -x0:dr:x0, pdf('Normal', -x0:dr:x0, 0, 1), 'r', ...
     R, 1./R.^alpha2, 'g');
grid on;
stmt = sprintf(['%s %dmin Return distribution %s -- %s'], ...
               strrep(company, '_', ' '), interval, start_day, ...
               end_day);
title(stmt);
legend('Empirical distribution', ...
       sprintf('power-law with xpnt %.4f', alpha1), ...
       'standard Gaussian', ...
       sprintf('power-law with xpnt %.4f', alpha2));
saveas(hdl, sprintf('pics/%s_%dmin_ret_%s-%s_lin.pdf', company, ...
                    interval, start_day, end_day));
% abs(prmt(2))), 'r-');
% plot(cen, empirical, '+');
