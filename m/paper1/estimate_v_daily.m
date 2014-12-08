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

R = NaN(T, p);
logvol = NaN(T, p);
c = 1;
for k = to_include'
    data = fetch(mysql, sprintf(...
        ['select closing, (high - low)/low from ' ...
         '%s_US order by day desc limit %d;'], ...
        strrep(symbols{k}, '.', '_'), T+1));
    R(:, c) = price2ret(flipud(cell2mat(data(:, 1))));
    A = flipud(cell2mat(data(2:end, 2)));
    A(A == 0) = min(A(A ~= 0));
    logvol(:, c) = log(A);
    logvol(:, c) = logvol(:, c) - mean(logvol(:, c));
    c = c + 1;
end
close(mysql);
