function ret = get_intra_ret(company, first_day, last_day, dt)
    days = get_recorded_days(company, first_day, last_day);

    num = ceil(60*8.5 - dt);

    ret = NaN(length(days) * num, 1);
    N = NaN(length(days), 1);
    n = 0;
    for d = 1:length(days)
        price = get_interpolated_prices(company, days{d}, 'minute', 1);
        if length(price) <= dt
            continue;
        end
        logprice = log(price);
        r = logprice(dt+1:end) - logprice(1:end-dt);
        N(d) = length(r);
        ret(n+1: n + N(d)) = r;
        n = n + N(d);
    end
    ret = ret(1:n);
    
