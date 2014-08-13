clear all
close all
T = 40;
horizon = 3;
day = '2014-07-18';
mysql = get_mysql();
stmt = sprintf(['select closing from DAX order by day desc limit ' ...
                '%d;'], T);
data = fetch(mysql, stmt);
close(mysql);

S = cell2mat(data(:, 1));
S = flipud(S);
r = price2ret(S);
