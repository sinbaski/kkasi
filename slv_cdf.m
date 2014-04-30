function y = slv_cdf(x, psi, sig)
x = reshape(x , length(x), 1);
da = 1.0e-2;
a = -10:da:10;
Ki = (ones(length(x), 1) * a) .* psi - x * exp(-sig .* a);
Ki = Ki ./ sqrt(2 * (1 - psi^2));
Ki = erfc(Ki);
Kr = exp(-a.^2 ./ 2) .* da;

y = Ki * Kr' ./ 2 ./ sqrt(2 * pi);
