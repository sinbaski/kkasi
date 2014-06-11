function mmt = johnson_su_moments34(args)
Gamma = args(1);
Delta = args(2);

w = exp(1 ./ Delta.^2);
Omega = Gamma ./ Delta;

m3 = sqrt(0.5 * w .* (w-1)) .* (...
     w .* (w + 2) .* sinh(3*Omega) + 3*sinh(Omega)) ./ ...
     (w .* cosh(2*Omega) + 1).^1.5;
m4 = (w.^2 .* (w.^4 + 2*w.^3 + 3*w.^2 - 3) .* cosh(4*Omega)...
      + 4 * w.^2 .* (w+2) .* cosh(2*Omega) + 3*(2*w + 1)) ./...
      (2 * (w .* cosh(2*Omega) + 1).^2);
mmt = [m3, m4];
