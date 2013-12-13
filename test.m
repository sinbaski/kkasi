clear all
filename = 'data/nasdaq.csv';
price = dlmread(filename, ',', [1, 4, 10807, 4]);
% filename = 'data/sp.csv';
% price = dlmread(filename, ',', [1, 4, 16084, 4]);
price = flipud(price); 
ret = price2ret(price);
ret = ret(abs(ret) < 5.0e-2);
r = (ret - mean(ret))/std(ret);

phi = fft(r);
l = length(r);

y = log(log(phi(2:end)) ./ log(phi(1:end-1)));
y = reshape(y, 1, l-1);
x = log([2:l] ./ [1:l-1]);
[alpha, b] = polyfit(x, y, 1);
