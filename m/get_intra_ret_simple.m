function [ret, sig] = get_intra_ret_simple(...
    company, first_day, last_day, dt, delta)

    days = get_recorded_days(company, first_day, last_day);

    num = 8.5*60/dt;
    ret = NaN(length(days) * num, 1);
    sig = NaN(length(days) * num, 1);

    % N = NaN(length(days), 1);
    % n = 0;
    last_closing = NaN;
    for d = 1:length(days)
        price = get_interpolated_prices(company, days{d}, 'minute', ...
                                                 dt);
        if isnan(last_closing)
            ret((d-1)*num+1 : d*num) = price2ret([price(1); price]);
        else
            ret((d-1)*num+1 : d*num) = price2ret([last_closing; price]);
        end
        last_closing = price(end);
        % if length(price) <= dt
        %     N(d) = 0;
        %     continue;
        % end
        % price = price(1:dt:length(price));
        % r = price2ret(price);
        % N(d) = length(r);
        % ret(n + 1: n + N(d)) = r;
        % n = n + N(d);
    end

    if isnan(delta) || isinf(abs(delta))
        return;
    end
    m = dt * 60 / delta;
    for d = 1:length(days)
        % A price observation every delta sec.
        price = get_interpolated_prices(...
            company, days{d}, 'second', delta);
        c = NaN;
        for k = 1:num
            % In each dt min, there are dt*2 intervals of 30s.
            r = price2ret(price((k-1)*m+1 : k*m));
            if norm(r)
                sig((d-1)*num+k) = norm(r);
                if (isnan(c))
                    c = k;
                end
            elseif k > 1
                sig((d-1)*num+k) = sig((d-1)*num+k-1);
            end
        end
        sig((d-1)*num+[1:c-1]) = sig((d-1)*num + c);
    end
    