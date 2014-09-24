% G: arguments to the Blue function, column vector, complex values
% v1: variance of the normal distribution, scalar
% ret: Values of the Blue function. Nx2 matrix. real values
function ret = LognormalBlue(G, v1, q1)
persistent Z;
persistent B;
persistent v q;

if (isempty(v) || v1 ~= v ||...
    isempty(q) || q1 ~= q ||...
    isempty(Z) || length(Z) ~= length(G) ...
    || sum(Z == G) ~= length(G))
    Z = G;
    v = v1;
    q = q1;
elseif Z == G
    ret = B;
    return;
end
    

L = 2000;
s = linspace(-6*sqrt(v), 6*sqrt(v), L);

M1 = 1./(1 - Z * exp(2*s) .* q);
M2 = exp(2.*s - s.^2/(2*v))'.*(s(2) - s(1));
F = (M1 * M2)./sqrt(2*pi*v) + 1./Z;
B = NaN(length(Z), 2);
B(:, 1) = real(F);
B(:, 2) = imag(F);
ret = B;
