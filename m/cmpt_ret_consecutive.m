function ret = cmpt_ret_consecutive(company, first_day, last_day, dt, width)

timefmt = 'yyyy-mm-dd HH:MM:SS';
hostname = '83.176.196.41';
%hostname = 'localhost';
jdbc = ['/nfs/users3/xiexiaol/lib/mysql-connector-java-5.1.26/'...
         'mysql-connector-java-5.1.26-bin.jar'];

if ~strcmp(hostname, 'localhost') && isempty(strfind(javaclasspath, jdbc))
    javaaddpath(jdbc);
end

mysql = database('avanza', 'sinbaski', 'q1w2e3r4',...
                 'com.mysql.jdbc.Driver', ...
                 sprintf('jdbc:mysql://%s:3306/avanza', hostname));
stmt = sprintf(['select distinct(date(tid)) from %s where date(tid)' ...
                '>= "%s" and date(tid) <= "%s" order by tid;'], company, ...
               first_day, last_day);
days = fetch(mysql, stmt);

ret = NaN(1000*length(days));
tot = 0;
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
    n1 = NaN;
    X = NaN(1000, 1);
    while m <= N - 1
        % N/102*dt: approx. # trades in dt
        stmt = sprintf(['select count(*)/timestampdiff(minute, min(tid), ' ...
                      'max(tid)) from %s where date(tid) ' ...
                      '= "%s";'], company, days{d});
        data = fetch(mysql, stmt);
        trade_per_min = cell2mat(data);
        l = dt/60 * trade_per_min;

        n1 = locateit(tid, m, l*dt, @intra_time_sort, dt);
        if isnan(n1) break; end
        if etime(tid(n1, :), tid(m, :)) > dt*(1 + width)
            m = n1;
            continue;
        end
        X(k) = log(price(n1)/price(m));
        k = k + 1;
        m = n1;
    end
    I = ~isnan(X);
    % Mu = mean(X(I, d));
    % Sigma = std(X(I, d));
    % X(I, d) = (X(I, d) - Mu) ./ Sigma;
    s = sum(I);
    ret(tot+1:tot+s) = X(I, d);
    tot = tot + s;
    fprintf('done with %s\n', char(days{d}));
end
ret = ret(~isnan(ret));
close(mysql);
