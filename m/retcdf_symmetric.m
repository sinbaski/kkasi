function p = retcdf_symmetric(x, s, evbar)
C = sqrt(2/pi)*exp(log(log(4))/2);
p = NaN(size(x));
I = isfinite(x);
x = x./evbar;

p(I) = (x(I)<0).*sqrt(log(4))./C.*exp(s^2./2)./sqrt(8.*pi).*(...
    x(I)./sqrt(log(4)).*erfc(...
        1./(s.*sqrt(2)) .* log(abs(x(I))./sqrt(log(4))) + s./sqrt(2))...
    + exp(-s^2./2).*erfc(1./(s.*sqrt(2)) .* log(abs(x(I))./sqrt(log(4)))))...
       +(x(I)>=0).*(1./2) + (x(I)>=0).*sqrt(log(4))./C.*exp(s^2./2)./sqrt(8.*pi).*(...
           x(I)./sqrt(log(4)).*erfc(...
               1./(s.*sqrt(2)) .* log(abs(x(I))./sqrt(log(4))) + s./sqrt(2))...
           + exp(-s^2./2).*erfc(-1./(s.*sqrt(2)) .* log(abs(x(I))./sqrt(log(4)))));


I = ~isfinite(x);
p(I) = (x(I) > 0);
