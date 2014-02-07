clear all
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


%% intraday returns
company = 'volvo_b';
first_day = '2013-10-10';
last_day = '2014-01-17';
dt = 30;
ret = get_intra_ret(company, first_day, last_day, dt);
r = ret - mean(ret);
sigma = std(r);

C = sqrt(2/pi)*exp(log(log(4))/2);
f = @(x, s, A) A*...
    exp(s^2/2)/sqrt(8*pi)/C .* ...
    erfc(1/(s*sqrt(2)) * log(abs(x*A)/sqrt(log(4))) + s/sqrt(2));

p = mle(r, 'pdf', f, 'start', [1, 1/std(r)], 'lowerbound', [1e-4, 1e-4]);
% [y, x] = ecdf(r);
% y1 = F(x, sig);
% plot(x, y, x, y1, 'r');

% x = linspace(min(r), max(r), length(r)/100);
% y = hist(r, x)./length(r)/(x(2)-x(1));
% plot(x, y, x, f(x, p(1), p(2)), 'r');


% plot(x(I), y(I), x(I), f(x(I), p(1), p(2)), 'r');
% plot(log(x(I)), log(y(I)), log(x(I)), log(f(x(I), p(1), p(2))),
% 'r');
[y, x] = ecdf(r);
y1 = retcdf(x, p(1), p(2));
I = x > 0 & x < sigma;
% I = x > 2*sigma & x < 3*sigma;
% I = logical(ones(1, length(x)));

plot(log(x(I)), log(1-y(I)), log(x(I)), log(1-y1(I)), 'r');
grid on
legend('ln(x) vs. ln(1 - ecdf(x))', 'ln(x) vs. ln(1 - F(x))', 'Location', 'Southwest');
title('Surviving probability distribution function');
xlabel('ln(x) with $x > 0$ and $x < \sigma_r$', 'Interpreter', 'latex', ...
       'Fontsize', 14);
ylabel('ln(1 - F(x)) and ln(1 - ecdf(x))');
