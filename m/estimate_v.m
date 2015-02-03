clear all
close all

company = 'nordea_bank';
dt = 30;
first_day = '2013-10-10';
last_day = '2014-08-21';
[ret, sig] = get_intra_ret_simple(company, '2013-10-10', '2014-08-21', ...
                                           dt, 60);
X = log(sig);
Z = ret ./ sig;

Z = Z .* exp(mean(X));
X = X - mean(X);


%% Try (high - low)/low as volatility proxy.
days = get_recorded_days(company, first_day, last_day);
mysql = get_mysql();
num = 8.5 * 60 / dt;
range_proxy = NaN(size(sig));
for d = 1:length(days)
    stmt = sprintf(['select min(price) as low , max(price) as high, ' ...
                    'floor(timestampdiff(MINUTE, "%s 09:00:00", tid)/%d) as ' ...
                    't from %s where date(tid) = "%s" group by t order by ' ...
                    't asc;'], days{d}, dt, company, days{d});
    data = cell2mat(fetch(mysql, stmt));
    A = NaN(num, 1);
    A(data(:, 3)+1) = (data(:, 2) - data(:, 1))./data(:, 1);
    c = NaN;
    for k = 2 : num
        if (isnan(A(k)))
            A(k) = A(k-1);
        elseif isnan(c)
            c = k;
        end
    end
    A(1:c-1) = A(c);
    A(A == 0) = mean(A ~= 0);
    range_proxy((d-1)*num+1 : d*num) = A;
end
close(mysql);
