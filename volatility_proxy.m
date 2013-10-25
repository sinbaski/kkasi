%clear all

% Volvo B 2007-06-01 -- 2008-12-31 down
% Volvo B 2009-01-01 -- 2010-12-31 up

filename = ['local_data/VOLV_B_SE0000115446-2003-09-12-2013-09-' ...
            '12.csv'];
filename = 'local_data/NDA-SEK-2012-01-16-2012-04-20.csv';
% opening  high  low  closing  average
price = dlmread(filename, ';', [1, 3, 2520, 7]);
price = flipud(price);
I = price(:, 2) ~= 0 & price(:, 3) ~= 0 & price(:, 4) ~= 0 & price(:, 5) ~= 0;

price = price(I, :);

p1 = price(1:912, :);
p2 = price(913:end, :);

vlt1(1:size(p1, 1)) = log(p1(:, 2) ./ p1(:, 3));
vlt2(1:size(p1, 1)) = log(p1(:, 2)./p1(:, 5));
vlt3(1:size(p1, 1)) = log(p1(:, 5)./p1(:, 3));
vlt4(1:size(p1, 1)) = (vlt1.^0.5 .* vlt2.^0.5);

vlt1(size(p1, 1)+1:size(p, 1)) = log(p2(:, 2) ./ p2(:, 3));
vlt2(size(p1, 1)+1:size(p, 1)) = log(p2(:, 2)./p2(:, 5));
vlt3(size(p1, 1)+1:size(p, 1)) = log(p2(:, 5)./p2(:, 3));
vlt4(size(p1, 1)+1:size(p, 1)) = (vlt1.^0.5 .* vlt2.^0.5);

%vlt5 = (vlt2 + vlt3) ./ 2;

[var(log(vlt1)), var(log(vlt2)), var(log(vlt3)), var(vlt4)]
