clear all
% close all
mysql = get_mysql();
data = cell2mat(fetch(mysql, ['select low, closing, high from DAX order by ' ...
                    'day desc limit 1000;']));
data = flipud(data);
close(mysql);

%% correlation with large returns
N = size(data, 1);
closing = data(1:end, 2);
ret = price2ret(closing);
% ranges = log(data(1:end, 3)) - log(data(1:end, 1));
Jup = NaN(N-59, 1);
Jdown = NaN(N-59, 1);
Short = NaN(N-59, 1);
c = 1;
d = 1;
e = 1;
A = [];
B = [];
% C = [];
for k = 60 : N-1
    if isempty(A)
        A = closing(k-59 : k);
        B = price2ret(A);
    else
        A = [A(2:end); closing(k)];
        B = [B(2:end); price2ret(closing(k-1:k))];
    end
    % C = price2ret(closing(k-9:9:k))/9;    
    if (sum(A <= closing(k)) <= 3 && ...
        ret(k) >= quantile(B, 0.9))
        Jup(c) = k+1;
        c = c + 1;
    elseif (ret(k) <= quantile(B, 0.1) &&...
            sum(A >= closing(k)) <= 3)
        Jdown(d) = k+1;
        d = d + 1;
    % elseif (sum(A > closing(k)) <= 3 &&...
    %         C > quantile(B, 0.75) &&...
    %         ret(k) < 0)
    %     Short(e) = k+1;
    %     e = e + 1;
    end
end
Jup = Jup(1 : c - 1);
Jdown = Jdown(1 : d - 1);
% Short = Short(1 : e-1);

plot(1:N, data(1:end, 2), Jup, data(Jup, 2), 'go',...
     Jdown, data(Jdown, 2), 'ro');
grid on

% N = size(data, 1);
% indicator = NaN(sum(I <= length(closing) - T3)-2, 1);
% target = NaN(sum(I <= length(closing) - T3)-2, 1);

% I = I(1:length(target)+2);
% for k = 3 : length(I)
%     indicator(k-2) = (ret(I(k)) -sum(ret(I(k-2:k-1))))/3;
%     target(k-2) = closing(I(k)+1+T3) - closing(I(k)+1);
% end

% plot(1:N, data(:, 2), I95, data(I95, 2), 'go', I5, data(I5, 2), 'ro');
% grid on


%% Correlation with ranges
% ret = price2ret(data(:, 2));
% ranges = log(data(:, 3)) - log(data(:, 1));
% closing = log(data(:, 2));
% T = 15;
% rhos = NaN(T, 1);
% levels = NaN(T, 1);
% for k = 1 : T
%     rhos(k) = corr(ranges(1:end-k), closing(k+1:end) - closing(k: ...
%                                                       end-1));
%     levels(k) = 2/sqrt(length(ranges)-k);
% end
% plot(1:T, rhos, 'o', 1:T, levels, 'r-', 1:T, -levels, 'r-');

%% correlation between returns of different time scales.
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
