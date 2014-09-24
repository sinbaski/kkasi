function den = LognormalDensity(lam, v, q, G)
a = erf(2*v);
b = 1 - (1/2) * exp(2*v) * erf(sqrt(8*v));
% h = (-a^3*b^3 - 3*a^2*b^4 + 15*a*b^5 - 10*b^6 - 27*b^5*q + ...
%      9*a*b^4*lam - 18*b^5*lam + sqrt(...
%          b^6*(-(a^2 + 2*a*b - 2*b*(b + 3*lam)).^3 + (a^3 + 3*a^2*b - ...
%                                                   3*a*b*(5*b + 3*lam) + ...
%                                                   b^2*(10*b + 27*q ...
%                                                   + 18*lam)).^2))).^(1/3);
       
% z = 2*b*(a+b) + b^2*(a^2 + 2*a*b-2*b*(b + 3*lam))./h + h;
% z = z ./ (6*b^2);

% z =
% 1/(6*b^2)*(2*b*(a+b)+(b^2*(a^2+2*a*b-2*b*(b+3*lam)))/(-a^3*b^3-3*a^2*b^4+15*a*b^5-10*b^6-27*b^5*q+9*a*b^4*lam-18*b^5*lam+sqrt(b^6*(-(a^2+2*a*b-2*b*(b+3*lam)).^3+(a^3+3*a^2*b-3*a*b*(5*b+3*lam)+b^2*(10*b+27*q+18*lam)).^2))).^(1/3)+(-a^3*b^3-3*a^2*b^4+15*a*b^5-10*b^6-27*b^5*q+9*a*b^4*lam-18*b^5*lam+sqrt(b^6*(-(a^2+2*a*b-2*b*(b+3*lam)).^3+(a^3+3*a^2*b-3*a*b*(5*b+3*lam)+b^2*(10*b+27*q+18*lam)).^2))).^(1/3));
% r = (q./(a-2*b*z) + 2*z - 1).^(-1/2);

options = optimset('Display', 'off');
f = @(arg) -2*b*arg^2 + (a + 2*b)*arg + b*q/(a - 2*b*arg) - b;
z = NaN(1, length(lam));

for k = 1:length(lam)
    z(k) = fsolve(@(arg) f(arg)-lam(k), 1, options);
end
r = (q./(a - 2*b*z) + 2*z - 1).^(-1/2);
den = (1/q/pi)*sqrt(r - r.^2 .* z.^2);
