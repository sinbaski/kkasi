% Study the auto-correlation
function NoneReturned = intra_autocorr(company, start_day, end_day, ...
                                       interval, round, get_data)

acf_datafile = sprintf('../local_data/%s_%dmin_fine_autocorr_%s-%s_%d_data.mat', company, ...
                       interval, start_day, end_day, round);

if get_data
    cmpt_intra_ret(company, start_day, end_day, interval);
end
if exist(acf_datafile, 'file') ~= 2
    cmpt_intra_acf(company, interval, start_day, end_day, round);
end

acf = [];
delta_t = [];
stmt = sprintf('ls ../local_data/%s_%dmin_fine_autocorr_%s-%s_*_data.mat', company, ...
               interval, start_day, end_day);
[status, output] = system(stmt);
files = strsplit(output);
for l = 1:length(files)
    load(files{l});
    acf = [acf, data.acf];
    delta_t = [delta_t, data.delta_t];
end

stmt = sprintf('ls ../data/%s_%dmin_ret_*.mat', company, interval);
[status, output] = system(stmt);
files = strsplit(output);

fmt = 'yyyy-mm-dd';
ds = datenum(start_day, fmt);
de = datenum(end_day, fmt);
%agrgt_ret = [];
done = 0;
% l = 1;
% while l <= length(files) && ~done
%     day = regexp(files(l), '[0-9]{4}-[0-9]{2}-[0-9]{2}', 'match');
%     day = datenum(day{1}, fmt);
%     if day >= ds && day <= de
%         load(files{l}, 'meta');
%         agrgt_ret = [agrgt_ret, meta.agrgt_ret];
%     end
%     if day == de
%         done = 1;
%     end
%     l = l + 1;
% end

T = 120 * 60;
%dt = T/800;
%cen = min(delta_t) : dt : T;
cen = 0:120;

% acf = [acf, agrgt_ret.^2];
% delta_t = [delta_t, zeros(1, length(agrgt_ret))];

% acf = (acf - mean(agrgt_ret)^2) ./ var(agrgt_ret, 1);

X = ones(1, length(cen)) * NaN;
%sigma = ones(length(cen), 1) * NaN;
for n = 1:length(cen)
    if n == 1
        I = delta_t == 0;
    else
        I = delta_t >= cen(n)*60 -  10& ...
            delta_t < cen(n)*60 + 10;
    end
    X(n) = mean(acf(I));
    %    sigma(n) = std(acf(I), 1);
end

% fprintf('max(acf) = %e, min(acf)=%e\n', max(X), min(X));
% save(sprintf('../local_data/%s_%dmin_fine_autocorr_%s-%s_std.txt', company, ...
%              interval, start_day, end_day), 'sigma', '-ascii');

hdl = figure;
plot(cen, X, 'b-');
grid on
title(sprintf(['%s %dmin ret. acf %s -- %s. %d acf val. %d bins. \n'...
               'min(acf)=%.4f  max(acf)=%.4f  acf(0)=%.4f)'], ...
              strrep(company, '_', ' '), interval, start_day, ...
              end_day, length(acf), length(cen), min(X), max(X), ...
              X(1)));

saveas(hdl, sprintf('../pics/%s_%dmin_fine_autocorr_%s-%s.pdf', company, ...
                    interval, start_day, end_day)); 

close(hdl);
