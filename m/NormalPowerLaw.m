function nloglf = NormalPowerLaw(param, x)
p = NaN(size(x));

xi1 = param(1);
u1 = sqrt(2*pi) * exp(xi1^2/2) * xi1 * normcdf(-xi1);
alpha1 = u1/(u1 - 1);

xi2 = param(2);
u2 = sqrt(2*pi) * exp(xi2^2/2) * xi2 * (1 - normcdf(xi2));
alpha2 = u2 / (u2 - 1);

I1 = x < -xi1;
I2 = x > xi2;
I3 = -xi1 <= x & x <= xi2;

p(I3) = normpdf(x(I3));

%% The left tail
a1 = (alpha1 / xi1)^alpha1 * exp(-xi1^2 / 2) / sqrt(2*pi);
mu1 = xi1 - alpha1/xi1;
if xi1 < mu1 || alpha1 < 1
    nloglf = Inf;
    return;
else
    p(I1) = a1 ./ abs(x(I1) + mu1).^alpha1;    
end



%% The right tail
a2 = (alpha2 / xi2)^alpha2 * exp(-xi2^2 / 2) / sqrt(2*pi);
mu2 = xi2 - alpha2/xi2;
if xi2 < mu2 || alpha2 < 1
    nloglf = Inf;
    return;
else
    p(I2) = a2 ./ (x(I2) - mu2).^alpha2;    
end

nloglf = -sum(log(p));
