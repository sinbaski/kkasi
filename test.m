clear all
close all

company = 'nordea_bank';
first_day = '2012-01-16';
last_day = '2012-04-20';

days = get_recorded_days('nordea_bank', '2012-01-16', ['2012-04-' ...
                    '20']);
ret15m = NaN(length(days) * 34, 1);
h15m = NaN(length(days) * 34, 1);

n = 1;
for d=1:length(days)
    price = get_interpolated_prices(company, days{d}, 'minute', 1);
    price = price(1:15:length(price));
    r = price2ret(price);
    N = length(r);
    ret15m(n: n + N - 1) = r;
    
    price = get_interpolated_prices(company, days{d}, 'second', 30);
    for k = 1:N
        r = price2ret(price((k-1) * 30 + 1 : k * 30));
        h15m(n+k-1) = sum(r.^2);
    end
    n = n + N;    
end
ret15m = ret15m(1:n-1);
h15m = h15m(1:n-1);

q = 5;
p = 1;
mdl0 = garch(p, q);
[mdl, covariance, loglikelihood, info] = estimate(mdl0, ret15m(q+1:end-1), ...
                                                  'E0', ret15m(1:q), ...
                                                  'V0', h15m(q-p+1:q));
% [mdl, covariance, loglikelihood, info] = estimate(mdl0, ret15m(q+1:end-1));
