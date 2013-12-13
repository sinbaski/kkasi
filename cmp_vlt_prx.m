clear all
close all
company = 'nordea_bank';
first_day = '2013-10-09';
last_day = '2013-12-12';
dt = 2;
% dt2 = 2;

% rv2 = cmpt_rlz_vlt(company, first_day, last_day, dt2);
% rv5 = cmpt_rlz_vlt(company, first_day, last_day, dt5);
filename = 'data/nordea.csv';
price = dlmread(filename, ',', [1, 2, 47, 5]);
price = flipud(price);
price = price([1:13, 15:size(price, 1)], :);
N = size(price, 1);

hl = log(price(:, 1) ./ price(:, 2));
volume = price(:, 4);
% volume = volume / mean(volume);

days = get_recorded_days(company, first_day, last_day);
ret = zeros(40, length(days));
for d = 1:length(days)
    S = get_interpolated_prices(company, days{d}, 'minute', 1);
    if length(S) <= dt
        continue;
    end
    S = S(1:dt:length(S));
    r = price2ret(S);
    ret(1:length(r), d) = r;
end
rv = sqrt(sum(ret.^2))';

P = polyfit(hl, rv, 1);
y = polyval(P, hl);

[hl1, I] = sort(hl);
rv1 = rv(I);

plot(hl1, rv1, 'x', hl1, polyval(P, hl1));

