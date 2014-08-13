%function rho = MarcenkoPasturPDF(lambda, q, var)
function rho = MarcenkoPasturPDF(lambda, prmt)
q = abs(real(prmt(1)));
var = abs(real(prmt(2)));
rho = NaN(size(lambda));

a = var * (1 - sqrt(q))^2;
b = var * (1 + sqrt(q))^2;

f = @(x) 1/(2*pi*q*var).*sqrt((b./x - 1) .* (1 - a./x));


fprintf('MarcenkoPasturPDF: q=%f, var=%f, a=%f, b=%f\n', ...
        q, var, a, b);

I = lambda > a & lambda < b;
rho(I) = f(lambda(I));

I1 = lambda < a;
I2 = lambda > b;

% rho(I1) = f(a + (b-a)/1000) .* exp(-(lambda(I1) - a).^2);
% rho(I2) = f(b - (b-a)/1000) .* exp(-(lambda(I2) - b).^2);
rho(I1 | I2) = 0;
%rho(~isreal(rho)) = 0;

