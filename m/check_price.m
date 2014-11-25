clear all

timefmt = 'yyyy-mm-dd HH:MM:SS';

start_day = '2012-03-16';
end_day = '2012-04-20';
company = 'Nordea_Bank';

mysql = database('avanza', 'sinbaski', 'q1w2e3r4',...
                 'com.mysql.jdbc.Driver', ...
                 'jdbc:mysql://localhost:3306/avanza');

stmt = sprintf(['select price from %s where date(tid)' ...
                '>= "%s" and date(tid) <= "%s" order by tid;'], company, ...
               start_day, end_day);

data = fetch(mysql, stmt);
close(mysql);
price = cell2mat(data(:, 1));
for k = 34100 : 34300
    if price(k) > 64; break; end
end
price(k) = [];
%price(34200) = [];
% m = mean(price);
% price = price(abs(price - m) / m < 0.2);

plot(1:length(price), price, 'b');
title(sprintf('%s %s -- %s paid prices', strrep(company, '_', ' '), ...
              start_day, end_day));
grid on;


