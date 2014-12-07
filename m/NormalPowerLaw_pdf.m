% function loglf = NormalPowerLaw_pdf(param, x)
function p = NormalPowerLaw_pdf(param, x)
p = NaN(size(x));

xi1 = param(1);
alpha1 = xi1 ./ normcdf(-xi1) .* exp(-xi1.^2 / 2) / sqrt(2*pi) + 1;
a1 = xi1^alpha1 * exp(-xi1^2 / 2) / sqrt(2*pi);

xi2 = param(2);
alpha2 = xi2 ./ (1 - normcdf(xi2)) .* exp(-xi2.^2 / 2) / sqrt(2*pi) + 1;
a2 = xi2^alpha2 * exp(-xi2^2 / 2) / sqrt(2*pi);

I1 = x < -xi1;
I2 = x > xi2;
I3 = -xi1 <= x & x <= xi2;

p(I3) = normpdf(x(I3));

%% The left tail

p(I1) = a1 ./ abs(x(I1)).^alpha1;    

%% The right tail

p(I2) = a2 ./ x(I2).^alpha2;    

%loglf = sum(log(p));