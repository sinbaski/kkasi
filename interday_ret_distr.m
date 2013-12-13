clear all
close all
% Volvo B 2007-06-01 -- 2008-12-31 down
% Volvo B 2009-01-01 -- 2010-12-31 up
% filename = 'data/volvo_b_20121203-20131202.csv';

%% daily returns
filename = 'data/nasdaq.csv';
price = dlmread(filename, ',', [1, 4, 10807, 4]);
% filename = 'data/sp.csv';
% price = dlmread(filename, ',', [1, 4, 16084, 4]);
price = flipud(price); 
ret = price2ret(price);
r = (ret - mean(ret))/std(ret);


%% intraday returns
% company = 'volvo_b';
% first_day = '2013-10-10';
% last_day = '2013-12-04';
% dt = 5;
% % interval for calculating realized volatility. in seconds
% delta = 10;
% [ret, v] = get_intra_ret_simple(...
%     company, first_day, last_day, dt, delta);
% r = (ret - mean(ret))/std(ret);

[p1, p2] = extract_xpnt(r, 2.56e-2/std(ret), 5.0e-2/std(ret));
% [p1, p2] = extract_xpnt(r, 0.045/std(ret));
fprintf('Left tail exponent: %.4f\nRight tail exponent: %.4f\n', ...
        p1, p2);

% figure;

% p = hist(r, x)/length(r)/(x(2) - x(1));
% plot(x, p);

% js = johnson_su_params([mean(r), var(r), skewness(r), ...
%                     kurtosis(r)]);
% pt = johnson_su_pdf(johnson_su_struct(js), x);
% plot(x, pt, 'g');


% param = mle(r, 'pdf',...
%             @(x, xi1, xi2)...
%             NormalPowerLaw_pdf([xi1, xi2], x),...
%             'start', sqrt([1.5, 1.42]), ...
%             'lowerbound', [1.0, 1.0]);
% powerlaw = NormalPowerLaw_cdf(param, x);

% figure;
% I1 = x < -param(1);
% plot(log(-x(I1)), log(p(I1)),...
%      log(-x(I1)), log(powerlaw(I1)));
% legend('empirical', 'power-law fit', ...
%        'Location', 'Southwest');
% title('S&P daily returns left tail');
% grid on;


% figure;
% I2 = x > param(2);
% plot(log(x(I2)), log(1-p(I2)),...
%      log(x(I2)), log(1-powerlaw(I2)));
% legend('empirical', 'power-law fit', ...
%        'Location', 'Northwest');
% title('S&P daily returns right tail');
% grid on



% [a1, a2] = extract_xpnt(ret);
% fprintf('a1 = %e, a2 = %e\n', a1, a2);


% cen = linspace(min(ret), max(ret), 80);
% dr = cen(2) - cen(1);
% empirical = double(hist(ret, cen)) ./ (length(ret) * dr);
% log_empirical = log(empirical);

% pd = fitdist(ret, 'Normal');

% subplot(3, 1, 1);
% plot(cen, log_empirical, 'b.', cen, log(pdf(pd, cen)), 'r');
% title('return distribution log-log scale');
% grid on;

% subplot(3, 1, 2);
% plot(cen, empirical, 'b.', cen, pdf(pd, cen), 'r');
% title('return distribution');
% grid on;

% subplot(3, 1, 3);
% ret = ret ./ stdret;

% acf = autocorr(ret, 15);
% plot(0:15, acf, 'bo-');

