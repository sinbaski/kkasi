% Z: arguments to the Blue function, column vector, complex values
% v: variance of the normal distribution, scalar
% D: Values of the Blue function's derivative. Nx2 matrix. real values
function D = LognormalBlueDer(Z, v, q)
L = 2000;
s = linspace(-6*sqrt(v), 6*sqrt(v), L);

M1 = 1./(1 - Z * exp(2*s) .* q).^2;
M2 = q.*exp(4.*s - s.^2/(2*v))'.*(s(2) - s(1));
F = (M1 * M2)./sqrt(2*pi*v) - 1./Z.^2;
D = NaN(2, 2, length(Z));
D(1, 1, :) = real(F);
D(1, 2, :) = -imag(F);
D(2, 1, :) = -D(1, 2, :);
D(2, 2, :) = D(1, 1, :);


