function y = fourierSumComplex(x, C, n, L)
k = [-n:0, 1:n]';
A = diag(C) * exp(-i * 2*pi * k * x ./ L);
y = sum(A);
