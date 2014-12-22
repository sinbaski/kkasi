clear all
mysql = get_mysql();
closing = cell2mat(fetch(mysql, ['select closing from DAX order by ' ...
                    'day desc limit 3000;']));
closing = log(flipud(closing));

ret = closing(2:end) - closing(1:end-1);

%% trend def. 1. irrelevant
% T1 = 10;
% T2 = 20;
% ind = NaN(length(ret)-T1+1, 1);
% for d = T1 : length(ret)
%     ind(d-T1+1) = sum(ret(d-T1+1:d) > 0)/T1 - 0.5;
% end

% ret20 = closing(T1+T2:end) - closing(T1:end-T2);

% ind = ind(1:length(ret20));


%% trend def. 2. irrelevant
% T1 = 10;
% T2 = 60;
% horizon = 20;
% ind = NaN(length(closing) - horizon - T2 +1, 1);
% for d = T2 : length(closing) - horizon
%     ind(d-T2+1) = mean(closing(d-T1+1:d)) - mean(closing(d-T2+1:d));
% end

% ret20 = closing(T2+horizon:end) - closing(T2:end-horizon);


%% 
ret20 = closing(21:end) - closing(1:end-20);
n = length(ret20);
corr(ret(2:1+n), ret20)
