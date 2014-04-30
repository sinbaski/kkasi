clear all
close all
% filename = './data/OMXS30.csv';
% filename = 'data/VOLV-B-1990-01-01-2014-01-09.csv';
% filename = 'data/ERIC-B-1990-01-02-2014-01-16.csv';
% data = dlmread(filename, ';', [1, 1, 6857, 3]);

% price = data(data(:, 1) > 0 & data(:, 3) > 0 &...
%             data(:, 1) ~= data(:, 2), 3);
% % price = data(data(:, 3) > 0, 3);
% price = flipud(price);
% N = size(price, 1);
% T = 15;
% % ret = price2ret(price(mod(N, T):T:end));
% % ret = ret - mean(ret);
% r = price2ret(price);
% r = r - mean(r);
% sigma = std(r);
% r = (r - mean(r))/std(r);
f1 = @(p) ...
     p(1) * p(2) * exp(p(2)^2 / 2 + p(3)) + p(4);

f2 = @(p) (4 * p(1)^2 * p(2)^2 + 1) * ...
     exp(2*p(3) + 2* p(2)^2) + 2*p(4)*p(2)*p(1)*exp(p(3) + p(2)^2/2)...
     + p(4)^2;

f3 = @(p) 9 * exp(3*p(3) + 4.5*p(2)^2) * p(2) ...
     * p(1) * (1 + 3*p(2)^2 * p(1)^2) + p(4)^3 ...
     + 3*exp(p(3) + p(2)^2/2)*p(4)^2*p(2)*p(1) ...
     + 3*exp(2*p(3) + 2*p(2)^2)*p(4)*(1 + 4*p(2)^2*p(1)^2);

f4 = @(p) exp(8*p(2)^2 + 4*p(3)) * ...
     (256 * p(1)^4 * p(2)^4 + 96 * p(1)^2 * p(2)^2 + 3)...
     + p(4)^4 + 4*exp(p(3) + p(2)^2/2)*p(4)^3*p(2)*p(1)...
     + 36*exp(3*p(3) + 9*p(2)^2/2)*p(4)*p(2)*p(1)*(1 + 3*p(2)^2*p(1)^2)...
     + 6*exp(2*p(3) + 2*p(2)^2)*p(4)^2*(1 + 4*p(2)^2*p(1)^2);

central = @(p) central_moments([f1(p), f2(p), f3(p), f4(p)]);

%% intraday returns
company = 'volvo_b';
first_day = '2013-10-10';
last_day = '2014-03-12';
dt = 30;
ret = get_intra_ret(company, first_day, last_day, dt);

noncentral = mean([ret, ret.^2, ret.^3, ret.^4]);

a = noncentral(3) / noncentral(2)^(3/2);
b = noncentral(4) /noncentral(2)^2;

c = a/b^(3/8);
f = @(x) 9*x .* (3*x.^2 + 1) ./ (4*x.^2 + 1).^(3/4) ./...
    (256 * x.^4 + 96*x.^2 + 3).^(3/8) - c;
x = -1e-1:1e-2:1e-1;
plot(x, f(x));
grid on

close all
y = fzero(f, 0);
sig = sqrt(2/3 * log((4*y ^2 + 1)^(3/2) / 9 / (3*y^2 + 1) * (a/y)));
psi = y / sig;
v = 1/2*log(noncentral(2)) - 1/2*log(4*y^2+1) - sig^2;

%% Estimate the parameters by matching the moments.
p0 = [psi, sig, v, 0];
p = fsolve(@(p) [f1(p); f2(p); f3(p); f4(p)], p0);

%% Verify the moments
% v = -2;
% sig = 0.4;
% psi = 0.6;
% mu = 0.1;

% p = [psi, sig, v, mu];
% X = mvnrnd([0, 0], [1, psi; psi, 1], 2e+5)';
% cov(X')
% r = exp(X(1, :) .* sig + v) .* X(2, :) + mu;

% [f1(p), mean(r)]
% [f2(p), mean(r.^2)]
% [f3(p), mean(r.^3)]
% [f4(p), mean(r.^4)]

%% Estimate the parameters by MLE
% q = mle(ret, 'pdf', @(x, p, s, v, m) exp(-v) .* slv_pdf(exp(-v).*(x-m), p, s), ...
%         'start', [p(1), p(2), p(3), p(4)], ...
%         'lowerbound', [-0.03, 0.3, -8, -0.01],...
%         'upperbound', [0.04, 0.6, log(max(abs(ret))), 0.01]);
% load('AsymmetricModelMLE.mat');
ops = statset('mlecustom');
ops.MaxFunEvals = 600;
ops.MaxIter = 300;

% AT LEAST ONE PARAMETER HAS TO BE FIXED, OTHERWISE MLE WON'T CONVERGE!!!
q = mle(ret, 'pdf', @(x, s, v, m) exp(-v) .* slv_pdf(exp(-v).*(x-m), p(1), s), ...
        'start', [p(2), p(3), p(4)], ...
        'lowerbound', [0.1, -8, -0.1],...
        'upperbound', [1, -6, 0.1], ...
        'options', ops);

%% Compare the CDF with empirical CDF using the estimated parameter values
q = [p(1), q];
M = central(q);
M(1) = f1(q);
M(2) = sqrt(M(2));
M(3) = M(3)/M(2)^3;
M(4) = M(4)/M(2)^4;
fprintf('MLE parameters:\n');
fprintf('%.4e & %.4e & %.4e & %.4e\n', q);
fprintf('Sample Moments:\n');
fprintf('%.4e & %.4e & %.4e & %.4e\n', mean(ret), std(ret), ...
        skewness(ret), kurtosis(ret));
fprintf('Theoretical Moments:\n');
fprintf('%.4e & %.4e & %.4e & %.4e\n', M);

[y, x] = ecdf(ret);
x1 = linspace(min(ret), max(ret), 500);
y1 = slv_cdf((x1 - q(4)).* exp(-q(3)), -0.0157, q(2));
% y2 = slv_cdf((x1 - p(4)).* exp(-p(3)), p(1), p(2));

subplot(1, 2, 1);
plot(log(-x(x<0)), log(y(x<0)), ...
     log(-x1(x1<0)), log(y1(x1<0)));
grid on
xlabel('ln(-x) where x < 0', 'Interpreter', 'tex');
ylabel('ln(P(r < x)) where x < 0', 'Interpreter', 'tex');
legend('empirical CDF', 'Model CDF', 'Location', 'Southwest');

subplot(1, 2, 2);
plot(log(x(x>0)), log(1-y(x>0)), log(x1(x1>0)), log(1-y1(x1>0)));
grid on
xlabel('ln(x) where x > 0', 'Interpreter', 'tex');
ylabel('ln(P(r > x)) where x > 0', 'Interpreter', 'tex');
legend('empirical compl. CDF', 'Model compl. CDF', 'Location', ...
       'Southwest');
