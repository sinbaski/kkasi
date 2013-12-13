%function rho = MarcenkoPasturPDF(lambda, q, var)
function rho = MarcenkoPasturPDF(prmt, lambda)
% a = N/T, i.e. # series / # observations 
q = abs(real(prmt(1)));
var = abs(real(prmt(2)));

a = var * (1 - sqrt(q))^2;
b = var * (1 + sqrt(q))^2;
fprintf('MarcenkoPasturPDF: q=%f, var=%f, a=%f, b=%f\n', ...
        q, var, a, b);

rho = 1/(2*pi*q*var) * ...
      sqrt((b./lambda - 1) .* (1 - a./ lambda));
rho(lambda < a | lambda > b) = 0;
%rho(~isreal(rho)) = 0;

