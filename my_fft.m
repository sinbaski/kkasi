function y = my_fft(x)
n = length(x);
coef = NaN(n, n);
for k = 1:n
    y(k) = 0;
    for j = 1:n
        y(k) = y(k) + exp(-i*2*pi*(k-1)*(j-1)/n) * x(j);
    end
end
