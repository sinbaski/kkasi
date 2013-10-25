function returns = get_intra_ret(company, first_day, last_day, ...
                                          interval, width, cmpt, ...
                                          agrgt)
if cmpt
    cmpt_intra_ret(company, first_day, last_day, interval, width);
end

stmt = sprintf('ls data/%s_%dsec_ret_*.mat', company, interval);
[status, output] = system(stmt);
files = strsplit(output);

fmt = 'yyyy-mm-dd';
ds = datenum(first_day, fmt);
de = datenum(last_day, fmt);
returns = [];
l = 1;
done = 0;
while l <= length(files) && ~done
    daystr = regexp(files(l), '[0-9]{4}-[0-9]{2}-[0-9]{2}', 'match');
    daynum = datenum(daystr{1}, fmt);
    if daynum >= ds && daynum <= de
        if agrgt
            load(files{l}, 'meta');
            returns = [returns; meta.agrgt_ret];
        else
            load(files{l}, 'ret');
            returns = [returns; ret];
        end
        if daynum == de
            done = 1;
        end
        l = l + 1;
    end
end
