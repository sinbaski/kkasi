% clear all
% func = @(x) x^2*(1 - log(x^2)) + log(2*pi)

% fzero(func, exp(0.5))
clear all

timefmt = 'yyyy-mm-dd HH:MM:SS';
interval = 20;
start_day = '2012-01-16';
end_day = '2012-03-15';
company = 'Nordea_Bank';

stmt = sprintf('ls data/%s_%dmin_ret_*.mat', company, interval);
[status, output] = system(stmt);
files = strsplit(output);
returns = [];
fmt = 'yyyy-mm-dd';
ds = datenum(start_day, fmt);
de = datenum(end_day, fmt);
l = 1;
done = 0;
data = struct('auto_corr', [], 'delta_t', [], ...
              'vol_auto_corr', []);
tot = 0;
while l <= length(files) && ~done
    day = regexp(files(l), '[0-9]{4}-[0-9]{2}-[0-9]{2}', 'match');
    day = datenum(day{1}, fmt);
    if day >= ds && day <= de
        load(files{l}, 'meta');
        tot = tot + length(meta.number);
    end
    if day == de
        done = 1;
    end
    l = l + 1;
end
fprintf('%d aggregated returns in total.\n', tot);
