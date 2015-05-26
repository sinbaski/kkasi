function [R, syms] = get_latest_data_from_index(...
    name, T, last_day, centering, operable)

mysql = get_mysql();

if strcmp(name, 'indices')
    stmt = sprintf('select name, tblname from %s;', name);
    data = fetch(mysql, stmt);
    tblname = data(:, 2);
    p = size(data, 1);
else
    if operable
        stmt = sprintf('select symbol from %s where operable;', name);
    else
        stmt = sprintf('select symbol from %s;', name);
    end
    symbols = fetch(mysql, stmt);
    tblname = symbols;
    
    p = size(symbols, 1);
    if strcmp(name, 'SP500_components')
        postfix = '_US';
    elseif strcmp(name, 'OMXS30_components')
        postfix = '_SE';
    else
        postfix = '';
    end
    for k = 1 : length(symbols)
        s = strrep(symbols{k}, '.', '_');
        s = strrep(s, '-', '_series_');
        s = strcat(s, postfix);
        tblname{k} = s;
    end
end

numRec = NaN(p, 1);
for k = 1 : p
    numRec(k) = cell2mat(fetch(mysql, sprintf(['select count(*) from %s'], tblname{k})));
end

to_include = find(numRec >= T+1);
if strcmp('indices', name)
    syms = data(to_include, 1);
else
    syms = symbols(to_include);
end

p = length(to_include);
R = NaN(T, p);
c = 1;
for k = to_include'
    closing = fetch(mysql, sprintf(['select closing from %s where day <= \"%s\" order ' ...
                        'by day desc limit %d;'], tblname{k}, last_day, T+1));
    R(:, c) = price2ret(cell2mat(flipud(closing)));
    if centering
        R(:, c) = R(:, c) - mean(R(:, c));
    end
    c = c + 1;
end
close(mysql);
