clear all
close all
% filename = 'data/sp.csv';
% price = dlmread(filename, ',', [1, 4, 16084, 4]);

% Volvo B 2007-06-01 -- 2008-12-31 down
% Volvo B 2009-01-01 -- 2010-12-31 up
% filename = 'data/volvo_b_20121203-20131202.csv';

%% daily returns
filename = 'data/OMXS30.csv';
price = dlmread(filename, ';', [1, 1, 6849, 6]);

price = price(price(:, 1) > 0 & price(:, 6) > 0, :);
price = price(price(:, 1) ~= price(:, 2), :);

flipud(price); 
ret = price2ret(price(:, 3));

v1 = log(price(:,1)./price(:, 2));
v1 = v1(2:end);

v2 = price(:, 6)/mean(price(:, 6));
v2 = v2(2:end);

corr([ret.^2, v1, v2])
% probplot(ret ./ sig);
