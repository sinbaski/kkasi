function code = cmpt_intra_ret(company, first_day, last_day, ...
                               interval, width)
timefmt = 'yyyy-mm-dd HH:MM:SS';
code = 0;
mysql = get_mysql();
days = get_recorded_days(company, first_day, last_day);
for d = 1:length(days)
    filename = sprintf(['data/%s_%dsec_ret_%s.mat'], ...
                       company, interval, days{d});
    stmt = sprintf(['select tid, price from %s where ' ...
                    'date(tid) = \"%s\" order by tid;'], company, ...
                   days{d});
    data = fetch(mysql, stmt);
    tid = datevec(data(:, 1), timefmt);
    price = cell2mat(data(:, 2));

    stmt = sprintf(['select count(*) from %s where ' ...
                    'date(tid) = \"%s\";'], company, days{d});
    data = fetch(mysql, stmt);
    tot = cell2mat(data);
    trade_per_min = tot / (60 * 8.5);
    l = interval/60 * trade_per_min;
    ret = [];
    t1 = interval*(1 - width);
    t2 = interval*(1 + width);
    mark1 = NaN;
    mark2 = NaN;
    N = length(price);
    meta = struct('time', [], 'agrgt_ret', [], 'number', [], ...
                  'volatility', []);
    for h = 1 : N - 1
        if isnan(mark1)
            mark1 = locateit(tid, h, l, @intra_time_sort, t1);
        end
        if isnan(mark1) break; end

        while (mark1 < N && ...
               etime(tid(mark1, :), tid(h, :)) < t1)
            mark1 = mark1 + 1;
        end
        if etime(tid(mark1, :), tid(h, :)) < t1 break; end
        if etime(tid(mark1, :), tid(h, :)) >= t2 continue; end
        if isnan(mark2) mark2 = mark1; end

        while (mark2 < N && ...
               etime(tid(mark2, :), tid(h, :)) < t2)
            mark2 = mark2 + 1;
        end
        if etime(tid(mark2, :), tid(h, :)) >= t2
            mark2 = mark2 - 1;
        end
        X = log(price(mark1:mark2) ./ price(h));
        ret = [ret; X];
        meta.agrgt_ret = [meta.agrgt_ret; mean(X)];
        meta.time = [meta.time; tid(h, :)];
        meta.number = [meta.number; mark2 - mark1 + 1];
    end
    save(filename, 'ret', 'meta');
    fprintf('Saved file %s\n', filename);
end
close(mysql);

