function price = get_interpolated_prices(company, day, unit, dt)

mysql = get_mysql();
stmt = sprintf(['select avg(price) as p, floor(timestampdiff(%s, ' ...
                '"%s 09:00:00", tid)/%d) as t from %s where date(tid) ' ...
                '= "%s" group by t;'], unit, day, dt, company, day);

data = fetch(mysql, stmt);
close(mysql);

p = cell2mat(data(:, 1));
t = cell2mat(data(:, 2)) + 1;

if strcmp(unit, 'minute') == 1
    N = 8.5*60/dt;
elseif strcmp(unit, 'second') == 1
    N = 8.5*3600/dt;
end

price = NaN(N, 1);
price(t) = p;

% interpolate. Potential out-of-range assignment in the previous
% line can result in zeros.

for k = t(1):t(end)
    if isnan(price(k)) || price(k) == 0
        price(k) = price(k-1);
    end
end
price = price(~isnan(price));

