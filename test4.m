clear all

rng('default');

%filename =
%'local_data/VOLV_B_SE0000115446-2003-09-12-2013-09-12.csv';

filename = 'local_data/NDA-SEK-2012-01-16-2012-04-20.csv';

% opening  high  low  closing  average
%price = dlmread(filename, ';', [1, 3, 2520, 7]);
price = dlmread(filename, ';', [1, 3, 68, 7]);

price = flipud(price);
I = price(:, 2) ~= 0 & price(:, 3) ~= 0 & price(:, 4) ~= 0 & price(:, 5) ~= 0;

price = price(I, :);

% p1 = price(1:912, :);
% p2 = price(913:end, :);
p1 = price;

vlt1 = log(p1(:, 2) ./ p1(:, 3));
vlt2 = log(p1(:, 2)./p1(:, 5));
vlt3 = log(p1(:, 5)./p1(:, 3));
vlt4 = (vlt1 .* vlt2 .* vlt3).^(1/3);

vlt = vlt1(2:end);
inno = price2ret(p1(:, 4));
inno = (inno - mean(inno)) ./ vlt;
%inno = (inno - mean(inno));

%probs = normcdf([-1.5, -0.5, 0.5, 1.5]);
q1 = [-2.1, -0.7, 0.7, 2.1];
probs = normcdf(q1);
q2 = quantile(inno, probs);
[R1, type1] = johnsrnd([q1; q2], 1e5, 1);
[F1, X1] = ecdf(R1);

% moments = {mean(inno), var(inno), skewness(inno), kurtosis(inno)};
% [R2, type2] = pearsrnd(moments{:}, 1e5, 1);
% [F2, X2] = ecdf(R2);

ecdf(inno);
hold on
stairs(X1, F1, 'r');
%stairs(X2, F2, 'g');
grid on
hold off;

% z = randn(1e6, 1);
% y = sinh(z);

% cen = linspace(min(y), max(y), 100);
% num = hist(y, cen);
% probden = double(num)/sum(num)/(cen(2)-cen(1));

% plot(cen, probden);

% [mean(y), var(y), skewness(y), kurtosis(y)]

