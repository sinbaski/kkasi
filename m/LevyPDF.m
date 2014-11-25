function y = LevyPDF(x, a, b, g, m)
k = linspace(-5/g, 5/g, 2000);
dk = k(2)-k(1);
x = reshape(x, length(x), 1);

A = exp(-i * x * k);
V = reshape(LevyFT(k, a, b, g, m), length(k), 1) .* dk;
y = A*V ./ (2*pi);




