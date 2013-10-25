clear all
company = 'Nordea_Bank';
interval = 30;
start_day = '2012-01-16';
end_day = '2012-03-15';

stmt = sprintf('ls data/%s_%dmin_ret_*.mat', company, interval);
[status, output] = system(stmt);
files = strsplit(output);
returns = [];
fmt = 'yyyy-mm-dd';
ds = datenum(start_day, fmt);
de = datenum(end_day, fmt);
l = 1;
done = 0;
data = struct('acf', [], 'delta_t', []);
while l <= length(files) && ~done
    day = regexp(files(l), '[0-9]{4}-[0-9]{2}-[0-9]{2}', 'match');
    day = datenum(day{1}, fmt);
    if day >= ds && day <= de
        load(files{l}, 'ret');
        returns = [returns; ret];
    end
    if day == de
        done = 1;
    end
    l = l + 1;
end
m1 = mean(returns)
m2 = var(returns)
