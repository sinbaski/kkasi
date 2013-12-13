function [ret, sig] = get_intra_ret_simple(...
    company, first_day, last_day, dt, delta)

    days = get_recorded_days(company, first_day, last_day);

    num = ceil(60/dt*8.5);

    ret = NaN(length(days) * num, 1);
    sig = NaN(length(days) * num, 1);
    N = NaN(length(days), 1);
    n = 1;
    for d = 1:length(days)
        price = get_interpolated_prices(company, days{d}, 'minute', 1);
        if length(price) <= dt
            continue;
        end
        price = price(1:dt:length(price));
        r = price2ret(price);
        N(d) = length(r);
        ret(n: n + N(d) - 1) = r;
        n = n + N(d);
    end
    ret = ret(1: n-1);

    if isnan(delta) || isinf(abs(delta))
        return;
    end
    n = 1;
    for d = 1:length(days)
        % A price observation every delta sec.
        price = get_interpolated_prices(...
            company, days{d}, 'second', delta);
        for k = 1:N(d)
            % In each dt min, there are dt*2 intervals of 30s.
            r = price2ret(price((k-1)*dt*60/delta+1 : k*dt*60/delta));
            sig(n+k-1) = norm(r);
        end
        n = n + N(d);
    end
    sig = sig(1:n-1);
