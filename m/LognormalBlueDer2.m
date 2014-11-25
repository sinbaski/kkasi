% Z1: arguments to the Blue function, column vector, complex values
% v1: variance of the normal distribution, scalar
% q1: Values of the Blue function. Nx2 matrix. real values
function ret = LognormalBlueDer2(Z1, v1, q1)
persistent Z;
persistent D;
persistent v q;

if (isempty(v) || v1 ~= v ||...
    isempty(q) || q1 ~= q ||...
    isempty(Z) || length(Z) ~= length(Z1) ...
    || sum(Z == Z1) ~= length(Z1))
    Z = Z1;
    v = v1;
    q = q1;
elseif Z == Z1
    ret = D;
    return;
end

L = 2000;
s = linspace(-6*sqrt(v), 6*sqrt(v), L);

M1 = 1./(1 - Z * exp(2*s) .* q).^3;
M2 = exp(6.*s - s.^2/(2*v))'.*(s(2) - s(1));
D = -2.*q^2.*(M1 * M2)./sqrt(2*pi*v) - 2./Z.^3;
ret = D;

     