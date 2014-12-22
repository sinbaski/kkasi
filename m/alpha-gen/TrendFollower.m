clear all
mysql = get_mysql();
closing = cell2mat(fetch(mysql, ['select closing from DAX order by ' ...
                    'day;']));
T1 = 60;
T2 = 10;
