function C = fourierExpandComplex(data, n, L)
N = length(data);
k = [-n:0, 1:n]';
data = reshape(data, 1, N);

A = exp(i * 2*pi * k * data ./ L);
C = A * ones(N, 1) ./ (N*L);

