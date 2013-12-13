clear all
filename = 'data/volvo_b.csv';
price = dlmread(filename, ',', [1, 4, 3616, 4]);
price = flipud(price); 

ret = price2ret(price);
r = (ret - mean(ret))/std(ret);

x = 1.0:0.02:2.4;
y = x;
n = length(x);

z = NaN(n, n);
for u = 1:n
    for v = 1:n
        z(u, v) = NormalPowerLaw_pdf([x(u), y(v)], r);
    end
end
[X, Y] = meshgrid(x);
figure;
surf(X, Y, z);


% u1 = sqrt(2*pi) * exp(xi1.^2/2) .* xi1 .* normcdf(-xi1);
% u1 = xi1 ./ normcdf(-xi1) .* exp(-xi1.^2 / 2) / sqrt(2*pi) + 1;
