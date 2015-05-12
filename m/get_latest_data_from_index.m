function R = get_latest_data_from_index(name, T, centering)
mysql = get_mysql();
stmt = sprintf('select symbol from %s_components;', name);
symbols = fetch(mysql, stmt);
p = size(symbols, 1);
if strcmp(name, 'SP500')
    postfix = '_US';
elseif strcmp(name, 'DAX')
    postfix = '_DE';
elseif strcmp(name, 'OMXS30')
    postfix = '_SE';
else
    postfix = '';
end

numRec = NaN(p, 1);
for k = 1 : p
    s = strrep(symbols{k}, '.', '_');
    s = strrep(s, '-', '_series_');
    numRec(k) = cell2mat(fetch(mysql, sprintf(['select count(*) from %s%s'], ...
                                              s, postfix)));
end

to_include = find(numRec >= T+1);
p = length(to_include);
R = NaN(T, p);
c = 1;
for k = to_include'
    s = strrep(symbols{k}, '.', '_');
    s = strrep(s, '-', '_series_');
    closing = fetch(mysql, sprintf(['select closing from %s%s order ' ...
                        'by day desc limit %d;'], s, postfix, T+1));
    R(:, c) = price2ret(cell2mat(flipud(closing)));
    if centering
        R(:, c) = R(:, c) - mean(R(:, c));
    end
    c = c + 1;
end
close(mysql);
