clear all
close all
mysql = get_mysql();
name='SP500';
stmt = sprintf('select symbol from %s_components;', name);
symbols = fetch(mysql, stmt);

p = size(symbols, 1);
numRec = NaN(p, 1);
day1 = '2000-01-02';
day2 = '2015-01-23';
for k = 1 : p
    data = fetch(mysql, sprintf(['select count(*) ' ...
                        'from %s_US where day between "%s" and "%s";'], ...
                                strrep(symbols{k}, '.', '_'), day1, ...
                                day2));
    numRec(k) = cell2mat(data(1, 1));
end
T = 3788;
to_include = find(numRec == T);

p = length(to_include);

R = NaN(T-1, p);
logvol = NaN(T - 1, p);
c = 1;
for k = to_include'
    data = fetch(mysql, sprintf(...
        ['select closing, (high - low)/low from ' ...
         '%s_US where day between "%s" and "%s" order by day;'], ...
        strrep(symbols{k}, '.', '_'), day1, day2));
    R(:, c) = price2ret(cell2mat(data(:, 1)));
    A = cell2mat(data(2:end, 2));
    % A(A == 0) = min(A(A ~= 0));
    logvol(:, c) = log(A);
    % logvol(:, c) = logvol(:, c) - mean(logvol(:, c));
    c = c + 1;
end
close(mysql);
