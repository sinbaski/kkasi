% Study the auto-correlation
clear all

timefmt = 'yyyy-mm-dd HH:MM:SS';
% in units of minute
interval = 30;
start_day = '2012-03-16';
end_day = '2012-04-20';
company = 'Nordea_Bank';

mysql = database('avanza', 'sinbaski', 'q1w2e3r4',...
                 'com.mysql.jdbc.Driver', ...
                 'jdbc:mysql://localhost:3306/avanza');
stmt = sprintf(['select distinct(date(tid)) from %s where date(tid)' ...
                '>= "%s" and date(tid) <= "%s" order by tid;'], company, ...
               start_day, end_day);
days = fetch(mysql, stmt);

N = 510 / interval;
m = 6;
acfs = ones(length(days), N) .* NaN;
acfs2 = ones(length(days), N) .* NaN;
dt = int32(interval/m);
R = ones(N, 1) .* NaN;
V = ones(N, 1) .* NaN;
C = zeros(N, N);
D = zeros(N, N);
for d = 1:length(days)
    start = fetch(mysql, sprintf(['select min(tid) from %s where date(tid) = ' ...
                        '"%s";'], company, days{d}));
    start = start{1};
    t = datenum(start, timefmt);
    for a = 1:m
        t = double((a-1)*dt) / (60*24) + t;
        stmt = sprintf(['select sum(price*volume)/sum(volume) ' ...
                        'from %s where date(tid) ' ...
                        '= "%s" and timestampdiff(SECOND, "%s", tid) ' ...
                        '>= %d group by floor(timestampdiff(SECOND, "%s", tid)/(60*%d)) ' ...
                        'order by tid;'], company, days{d}, start, ...
                       (a-1)*dt*60, datestr(t, timefmt), interval);
        data = fetch(mysql, stmt);
        price = cell2mat(data(:, 1));
        ret = price2ret(price);
        vol = abs(ret);
        
        R = (ret - mean(ret)) ./ std(ret);
        V = (vol - mean(vol)) ./ std(vol);
        C = C + bsxfun(@times, R, R');
        D = D + bsxfun(@times, V, V');
        % n = min(N, length(ret)-1);
        % acfs((d-1)*m+a, 1:n) = autocorr(ret, n-1)';
        % acfs2((d-1)*m+a, 1:n) = autocorr(abs(ret), n-1)';
    end
    fprintf('Done with %s\n', days{d});
end
close(mysql);

X = C ./ (length(days) * m);
acfs = X(1, :);

X = D./N;
acfs2 = X(1, :);

hdl = figure;
plot(0:N-1, mean(acfs), 'bx-', 0:N-1, mean(acfs2), 'r+-');
stmt = sprintf(['%s %dmin Return/volatility autocorrelation %s -- %s'], ...
               strrep(company, '_', ' '), interval, start_day, ...
               end_day);
title(stmt);
legend('<r_i r_j>', '<|r_i| |r_j|>');
grid on
saveas(hdl, sprintf('pics/%s_%dmin_ret_%s-%s_autocorr.pdf', company, ...
                        interval, start_day, end_day));
