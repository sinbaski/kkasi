function y = InverseFourerTransform(x, k, fk)
k = reshape(k, 1, length(k));
x = reshape(x, length(x), 1);
fk = reshape(fk, length(fk), 1);

dk = k(2) - k(1);

A = exp(-i * x * k);
y = A * fk .* dk / (2*pi);

