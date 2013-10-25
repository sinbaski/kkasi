clear all
v = (exp(1)^2 - 1)/2;
x = -4.5:1e-3:4.5;
y_n = pdf('Normal', x, 0, v^0.5);
y_su = (2*pi*(1 + x.^2)).^(-1/2) .* exp(-asinh(x).^2 ./ 2);
plot(x, y_n, 'b', x, y_su, 'r');

