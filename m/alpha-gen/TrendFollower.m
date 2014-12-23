clear all
close all
mysql = get_mysql();
data = cell2mat(fetch(mysql, ['select low, closing, high from DAX order by ' ...
                    'day desc limit 250;']));
data = flipud(data);
close(mysql);

% ret = price2ret(data(:, 2));
ranges = log(data(:, 3)) - log(data(:, 1));
closing = log(data(:, 2));
T = 15;
rhos = NaN(T, 1);
levels = NaN(T, 1);
for k = 1 : T
    rhos(k) = corr(ranges(1:end-k), closing(k+1:end) - closing(k: ...
                                                      end-1));
    levels(k) = 2/sqrt(length(ranges)-k);
end
plot(1:T, rhos, 'o', 1:T, levels, 'r-', 1:T, -levels, 'r-');




% closing = log(flipud(closing));
% T = 60;
% horizons = [1, 5, 10, 15, 20, 25];
% ret1 = closing(2:end) - closing(1:end-1);
% rhos = NaN(T-1, 1);
% for m = 1 : length(horizons)
%     h = horizons(m);
%     retH = closing(T+h : end) - closing(T : end-h);
%     N = length(retH);

%     for k = 1 : T-1
%         rhos(k) = corr(ret1(T-k : T-k+N-1), retH);
%     end
%     subplot(3, 2, m);
%     plot(1:T-1, rhos, 'o', ...
%          1:T-1, 2*ones(1, T-1)./sqrt(N), ...
%          1:T-1, -2*ones(1, T-1)./sqrt(N));
%     title(sprintf('%d day horizon', h));
% end
