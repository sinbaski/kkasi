clear all
close all
% company = 'nordea_bank';
% first_day = '2012-01-16';
% last_day = '2012-04-30';
% dt5 = 5;
% dt2 = 2;

% rv2 = cmpt_rlz_vlt(company, first_day, last_day, dt2);
% rv5 = cmpt_rlz_vlt(company, first_day, last_day, dt5);
load('data/nordea_rv.mat', 'rv2', 'rv5');
filename = 'data/nordea_bank.csv';
price = dlmread(filename, ',', [369, 2, 436, 4]);
volume = dlmread(filename, ',', [369, 5, 436, 5]);
N = size(price, 1);
I = ones(N, 1);
I(13) = 0;
I(19) = 0;
I(37) = 0;
I(53) = 0;
price = price(logical(I), :);
volume = volume(logical(I));

price = flipud(price);
volume = flipud(volume);

hl = log(price(:, 2) ./ price(:, 3));
% x = (rv2 - mean(rv2)) .* (volume - mean(volume)) ./ (std(rv2) *
% std(volume));
y = ((rv2 - mean(rv2)) .* (hl - mean(hl))) ./ (std(rv2) * std(hl));

% x = (x - mean(x))./std(x);

% y = hl./rv2;
% y = (y - mean(y)) ./ std(y);
