clear all
close all
timefmt = 'yyyy-mm-dd HH:MM:SS';
% in units of minute
dt = 45;
start_day = '2012-01-16';
end_day = '2012-04-20';
company = 'nordea_bank';

[ret, sig] = get_intra_ret_simple(company, '2012-01-16', '2012-04-20', dt);
r = (ret - mean(ret))./sig;

q1 = [-2.1, -0.7, 0.7, 2.1];
probs = cdf('Normal', q1, 0, 1);
q2 = quantile(r, probs);
[R1, type1] = johnsrnd([q1; q2], 1e5, 1);
[F1, X1] = ecdf(R1);

ecdf(r);
hold on
stairs(X1, F1, 'r');
%stairs(X2, F2, 'g');
grid on
hold off;

sk = 0.1264;
kts = 3.1940;

f1 = @(x) - sqrt(0.5 * x(1) .* (x(1)-1)) .* (...
     x(1) .* (x(1) + 2) .* sinh(3*x(2)) + 3*sinh(x(2))) ./ ...
     (x(1) .* cosh(2*x(2)) + 1).^1.5;

f2 = @(x) (x(1).^2 .* (x(1).^4 + 2*x(1).^3 + 3*x(1).^2 - 3) .* cosh(4*x(2))...
      + 4 * x(1).^2 .* (x(1)+2) .* cosh(2*x(2)) + 3*(2*x(1) + 1)) ./...
      (2 * (x(1) .* cosh(2*x(2)) + 1).^2);

f = @(x) [f1(x) - sk; f2(x) - kts];

g = -0.7373;
d = 4.787;

X = fsolve(f, [exp(1/d^2); g/d]);

d = log(X(1)).^(-1/2);
g = X(2) .* d;

