function NoneReturned = cmpt_acf2(company, start_day, end_day, ...
                                  interval)
timefmt = 'yyyy-mm-dd HH:MM:SS';
hostname = '83.176.196.41';

mysql = database('avanza', 'sinbaski', 'q1w2e3r4',...
                 'com.mysql.jdbc.Driver', ...
                 sprintf('jdbc:mysql://%s:3306/avanza', hostname));
stmt = sprintf(['select distinct(date(tid)) from %s where date(tid)' ...
                '>= "%s" and date(tid) <= "%s" order by tid;'], company, ...
               start_day, end_day);
days = fetch(mysql, stmt);

ret = NaN(200, length(days));
for d = 1:length(days)
    stmt = sprintf(['select tid, price from %s where ' ...
                    'date(tid) = \"%s\" order by tid;'], company, ...
                   days{d});
    data = fetch(mysql, stmt);
    N = size(data, 1);
    tid = datevec(data(:, 1), timefmt);
    price = cell2mat(data(:, 2));
    m = 1;
    k = 1;
    while m <= N - 1
        % N/102*interval: approx. # trades in interval
        n = locateit(tid, m, N/30600*interval, @intra_time_sort, ...
                     interval);
        if isnan(n)
            break;
        end
        if etime(tid(n, :), tid(m, :)) > interval*1.03
            m = n;
            continue;
        end
        ret(k, d) = log(price(n)/price(m));
        k = k + 1;
        m = n;
    end
    I = ~isnan(ret(:, d));
    Mu = mean(ret(I, d));
    Sigma = std(ret(I, d));
    ret(I, d) = (ret(I, d) - Mu) ./ Sigma;
    fprintf('done with %s\n', char(days{d}));
end
close(mysql);

% 2 hours
numLag = 5;

Acf = NaN(numLag+1, size(ret, 2));
for d = 1:size(ret, 2)
    I = ~isnan(ret(:, d));
    Acf(:, d) = autocorr(ret(I, d), numLag);
end
Y = mean(Acf, 2);

hdl = figure;
plot((0:numLag)*interval, Y, 'bo-');
grid on;
xlabel('time lag (second)');
title(sprintf('%s %dsec ret. ACF %s -- %s', strrep(company, '_', ' '), ...
              interval, start_day, end_day));
saveas(hdl, sprintf('../pics/%s_%dsec_autocorr_%s-%s.pdf', company, ...
                    interval, start_day, end_day)); 

close(hdl);
