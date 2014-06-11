clear all
mysql = get_mysql();
stmt = sprintf(['select price from nordea_bank where date(tid) ' ...
                'between "2012-01-16" and "2012-04-20" order by tid;' ...
                '']);

data = fetch(mysql, stmt);
prices = cell2mat(data(:, 1));
close(mysql);
