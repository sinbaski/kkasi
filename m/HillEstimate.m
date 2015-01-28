clear all
close all
mysql = get_mysql();
name='SP500';
stmt = sprintf('select symbol from %s_components;', name);
symbols = fetch(mysql, stmt);

p = size(symbols, 1);
T = 3000;
numRec = NaN(p, 1);

for k = 1 : p
    numRec(k) = cell2mat(fetch(mysql, sprintf(['select count(*) from %s_US'], ...
                                              strrep(symbols{k}, '.', '_'))));
end

to_include = find(numRec >= T+1);
p = length(to_include);

tail_index = NaN(p, 2);
c = 1;
names = {};
for k = to_include'
    closing = fetch(mysql, sprintf(['select closing from %s_US order by day desc ' ...
                        'limit %d;'], strrep(symbols{k}, '.', '_'), ...
          T+1));
    R = price2ret(cell2mat(flipud(closing)));
    
    qlower = quantile(R, 0.03);
    Rlower = R(R < qlower);
    tail_index(c, 1) = 1/mean(log(Rlower./qlower));
    
    qupper = quantile(R, 0.97);
    Rupper = R(R > qupper);
    tail_index(c, 2) = 1/mean(log(Rupper./qupper));
    names{c} = symbols{k};
    c = c + 1;
end
close(mysql);
plot(tail_index(:, 2), tail_index(:, 1), '.');
hold on
X = linspace(0, 4, 400);
plot(X, X, 'r', 'Linewidth', 2);
xlabel('Upper tail index');
ylabel('lower tail index');
grid on
