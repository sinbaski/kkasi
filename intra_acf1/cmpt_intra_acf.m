function code = cmpt_intra_acf(company, interval, start_day, end_day, ...
                               round)
acf_datafile = sprintf('../local_data/%s_%dmin_fine_autocorr_%s-%s_%d_data.mat', company, ...
                       interval, start_day, end_day, round);

stmt = sprintf('ls ../data/%s_%dmin_ret_*.mat', company, interval);
[status, output] = system(stmt);
files = strsplit(output);

fmt = 'yyyy-mm-dd';
ds = datenum(start_day, fmt);
de = datenum(end_day, fmt);
data = struct('acf', [], 'delta_t', []);
T = 120;
done = 0;

variance = NaN * ones(1, length(files));
mu = NaN * ones(1, length(files));

if round == 0
    l = 1;
    while l <= length(files) && ~done
        daystr = regexp(files(l), '[0-9]{4}-[0-9]{2}-[0-9]{2}', 'match');
        daynum = datenum(daystr{1}, fmt);
        if daynum >= ds && daynum <= de
            load(files{l});
            mu = mean(meta.agrgt_ret);
            variance = var(meta.agrgt_ret, 1);
            data.acf = (meta.agrgt_ret - mu).^2 ./ variance;
            data.delta_t = zeros(1, length(meta.agrgt_ret));
            save(sprintf('../local_data/%s_mean_n_var_%s.mat', company, ...
                         char(daystr{1})), 'mu', 'variance');
        end
        l = l + 1;
    end
end

l = 1;
while l <= length(files) && ~done
    daystr = regexp(files(l), '[0-9]{4}-[0-9]{2}-[0-9]{2}', 'match');
    daynum = datenum(daystr{1}, fmt);
    if daynum >= ds && daynum <= de
        load(files{l});
        load(sprintf('../local_data/%s_mean_n_var_%s.mat', company, ...
                     char(daystr{1})));
        N = length(meta.agrgt_ret);
        k = 1;
        for m = 1:N-1
            for n = m+1:N
                data.acf(k) = (meta.agrgt_ret(m) - mu)*(meta.agrgt_ret(n) - mu)...
                    / variance;
                data.delta_t(k) = etime(meta.time(n, :), meta.time(m, ...
                                                                  :));
                k = k + 1;
            end
        end
        
        % q = int32(sqrt(N)*sqrt(2));
        % mark1 = round*q + 1;
        % mark2 = (round+1)*q;
        % if mark2 > N
        %     fprintf('passed %s\n', files{l});
        %     l = l + 1;
        %     continue;
        % end
        % U = meta.agrgt_ret(mark1 : mark2) - mu;
        % tau = etime(meta.time(end, :), meta.time(mark2, :));
        % if tau == 0
        %     s = N;
        % else
        %     distance = (N - mark2)*T*60/tau;
        %     s = locateit(meta.time, mark2, distance, @intra_time_sort, T*60);
        %     if isnan(s) s = N; end
        % end
        % for t = mark2 + 1 : s
        %     V = meta.agrgt_ret(t - q + 1 : t) - mu;
        %     data.acf = [data.acf, U.*V ./ variance];
        %     delta_t = abs(etime(meta.time(t-q+1 : t, :), ...
        %                         meta.time(mark1 : mark2, :)));
        %     data.delta_t = [data.delta_t, delta_t'];
        % end
        fprintf('processed %s\n', files{l});
    % else
    %     fprintf('ignored %s\n', files{l});
    end
    if daynum == de
        done = 1;
    end
    l = l + 1;
end
save(acf_datafile, 'data');
