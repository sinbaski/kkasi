clear all
filename = './data/OMXS30.csv';
data = dlmread(filename, ';', [1, 1, 6849, 6]);

% data = data(data(:, 1) > 0 && data(:, 6) > 0 &&...
%             data(:, 1) ~= data(:, 2), :);
price = data(data(:, 3) > 0, 3);
price = flipud(price);
N = size(price, 1);

T = 15;
ret = price2ret(price(mod(N, T):T:end));
mu = mean(ret);
r = price2ret(price(mod(N, T):1:end));
r = r - mean(r);
r = reshape(r, T, length(ret));
sig = sqrt(sum(r.^2))';

logv = log(sig);
L = length(logv);
rho = reshape(autocorr(logv, L-1), 1, L);

a = nan(1, L);
omg = nan(1, L);
for n = 0:L-1
    a(n+1) = 2*(rho./sqrt(L))*(cos(2*pi*n*[0:L-1]/L)'./sqrt(L));
end
omg = 2*pi*[0:L-1]/L;
I = omg < 1;
plot(omg(I), log(a(I)));
% plot(log(omg(I)), log(a(I)), 'x');

