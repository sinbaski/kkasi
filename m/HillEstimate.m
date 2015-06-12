clear all
close all
mysql = get_mysql();
name='SP500';
stmt = sprintf('select symbol from %s_components;', name);
symbols = fetch(mysql, stmt);

p = size(symbols, 1);
T = 1600;
numRec = NaN(p, 1);

for k = 1 : p
    name = strrep(symbols{k}, '.', '_');
    name = strrep(name, '-', '_series_');
    numRec(k) = cell2mat(fetch(mysql, sprintf(['select count(*) ' ...
                        'from %s_US'], name)));
end

to_include = find(numRec >= T+1);
p = length(to_include);

tail_index = cell(p, 3);
c = 1;
for k = to_include'
    name = strrep(symbols{k}, '.', '_');
    name = strrep(name, '-', '_series_');

    closing = fetch(mysql, sprintf(['select closing from %s_US order ' ...
                        'by day desc limit %d;'], name, T+1));
    R = price2ret(cell2mat(flipud(closing)));
    
    tail_index{c, 1} = symbols{k};
    
    qlower = quantile(R, 0.03);
    Rlower = R(R < qlower);
    tail_index{c, 2} = 1/mean(log(Rlower./qlower));
    
    qupper = quantile(R, 0.97);
    Rupper = R(R > qupper);
    tail_index{c, 3} = 1/mean(log(Rupper./qupper));

    c = c + 1;
end
close(mysql);

save('SP500TailIndices.mat', 'tail_index');
plot(cell2mat(tail_index(:, 3)), cell2mat(tail_index(:, 2)), '.');
hold on
X = linspace(0, 4, 400);
plot(X, X, 'r', 'Linewidth', 2);
xlabel('Upper tail index');
ylabel('lower tail index');
title('Tail indices of 441 S&P500 stocks');
grid on
