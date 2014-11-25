clear all
close all
mysql = get_mysql();
stmt = sprintf('select symbol from DAX_components;');
data = fetch(mysql, stmt);

p = 29;
T = 3000;

R = NaN(p, T);
for k = 1 : p
    closing = fetch(mysql, sprintf(['select closing from %s order by day desc ' ...
                        'limit %d;'], strrep(data{k}, '.', '_'), ...
          T));
    R(k, :) = cell2mat(flipud(closing))';
end
close(mysql);

