function [a, b] = fourierExpand(x, y, dx, x0, x1, n)
L = (x1 - x0);
r = [1:n]';

x = reshape(x, 1, length(x));
y = reshape(y, 1, length(y));
A = 2*pi*r * x ./ L;
a = cos(A) * (y' .* dx) .* 2 ./ L;
b = sin(A) * (y' .* dx) .* 2 ./ L;
