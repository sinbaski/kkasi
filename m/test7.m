clear all

filename = './data/OMXS30.csv';
data = dlmread(filename, ';', [1, 1, 6849, 6]);

% data = data(data(:, 1) > 0 && data(:, 6) > 0 &&...
%             data(:, 1) ~= data(:, 2), :);
price = data(data(:, 3) > 0, 3);
price = flipud(price);
N = size(price, 1);

T = 25;
ret = price2ret(price(mod(N, T):T:end));
mu = mean(ret);
r = price2ret(price(mod(N, T):1:end));
r = r - mean(r);
r = reshape(r, T, length(ret));
sig = sqrt(sum(r.^2))';

logv = log(sig);

w = ts_difference(logv, 15, 10);

coef = zeros(45, 45);
confidence = zeros(45, 45);
for T1 = 1:45
    for T2 = 1:45
        r1 = log(price(T1+1:end-T2)./price(1:end-T2-T1));
        r1 = r1-mean(r1);
        r2 = log(price(T1+1+T2:end)./price(T1+1:end-T2));
        r2 = r2 - mean(r2);
        coef(T1, T2) = corr(r1, r2);
        if abs(coef(T1, T2)) < 2/sqrt(N-T1-T2)
            coef(T1, T2) = 0;
        end
    end
end
[x, y] = meshgrid([1:45], [1:45]);
figure;
surf(x, y, coef);
xlabel('T2');
ylabel('T1');
% hold on
% surf(x, y, confidence, 'Facecolor', 'black', 'Edgecolor', 'black');
% hold off

