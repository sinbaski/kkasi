function ret = get_intra_ret_simple(company, first_day, last_day, dt)

    days = get_recorded_days(company, first_day, last_day);

    ret = NaN(1e6, 1);
    m = 1;
    for d = 1:length(days)
        price = get_interpolated_prices(company, days{d}, 'minute', 1);
        N = length(price);
        if N <= dt
            continue;
        end
        ret(m : m + N - dt - 1) = ...
            log(price(1 + dt : end) ./ price(1:end - dt));
        m = m + (N - dt);
    end
    ret = ret(1:m-1);

