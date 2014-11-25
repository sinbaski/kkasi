clear all
close all

company = 'volvo_b';
first_day = {
    '2014-02-14'
    '2014-01-31'
    '2014-01-17'
    '2014-01-03'
    '2013-12-20'
    '2013-12-06'
    '2013-11-22'
    '2013-11-08'
            };

last_day = {
    '2014-07-13'
    '2014-06-29'
    '2014-06-15'
    '2014-06-01'
    '2014-05-06'
    '2014-04-22'
    '2014-04-08'
    '2014-03-25'
           };
dt = 15;
delta = 30;
lags = 500;
acf = NaN(size(first_day, 1), 1000, 20);
ret = NaN(size(first_day, 1), 5000);
texts = {};
hold on
% for k = 1 : size(first_day, 1)
%     [r, v] = get_intra_ret_simple(...
%         company, first_day{k}, last_day{k}, dt, delta ...
%         );
%     ret(k, 1:length(r)) = r;
%     acf(k, 1:lags+1) = autocorr(r, lags);
%     plot(1:lags, acf(k, 2:lags+1));
%     texts{k} = sprintf('%s - %s', first_day{k}, last_day{k});
%     fprintf('Done with No. %d.\n', k);
% end

lags = 1;
p = 1;
acf = NaN(size(first_day, 1), 2, 81);
for nu = -4:0.1:4
    for k = 1 : size(first_day, 1)
        r = ret(k, :);
        r = r(~isnan(r));
        acf(k, 1:lags+1, p) = autocorr(r.^nu, lags);
    end
    p = p +1;
end
nu = -4:0.1:4;
plot((1:lags)', (acf(1:2:7, 2:lags+1)'));
ylabel('Autocorrelation', 'Fontsize', 14);
xlabel('#lags', 'Fontsize', 14);
legend(texts{1:2:7});

