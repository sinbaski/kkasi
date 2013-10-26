mysql = get_mysql();
data = fetch(mysql, ['select price from nordea_bank where date(tid) '...
                    'between "2012-01-16" and "2012-04-20" '...
                    'order by tid;']);
P = cell2mat(data(:, 1));
plot(1:length(P), P);
grid on
xlabel('# of trades');
ylabel('price (SEK)');
