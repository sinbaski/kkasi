function y = ecf(data, k)
y = NaN(size(k));
for n = 1:length(k)
    y(n) = sum(exp(i * k(n) .* data));
end
y = y ./ length(data);
    