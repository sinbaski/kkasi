function vlt = cmpt_rlz_vlt(company, first_day, last_day, dt)
timefmt = 'yyyy-mm-dd HH:MM:SS';
days = get_recorded_days(company, first_day, last_day);

vlt = NaN(length(days), 1);
for d = 1:length(days)
    price = get_interpolated_prices(company, days{d}, 'minute', 1);
    N = length(price);
    X = zeros(N, dt);
    Y = zeros(1, dt);
    for n = 1 : size(X, 2)
        l = length(n:dt:N) - 1;
        X(1:l, n) = price2ret(price(n:dt:N));
        Y(n) = sqrt(sum(X(1:l, n).^2));
    end
    vlt(d) = mean(Y);
    fprintf('rv%d for %s is %.4f\n', dt, days{d}, vlt(d));
end

