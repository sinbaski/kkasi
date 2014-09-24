% v: scalar
% q: scalar
% a: scalar or column vector
% return: scalar or column vector
function [f, fder] = LognormalEndsFun1(v, q, a)
L = 4000;
s = linspace(-5*sqrt(v), 5*sqrt(v), L);
M1 = exp(-s.^2./2./v) .* (s(2) - s(1));
M2 = repmat(exp(-2.*s), length(a), 1) - q.*repmat(a, 1, L);

M3 = 1./M2.^2;
f = M3 * M1' .* (q / sqrt(2*pi*v)); 

M4 = 1./M2.^3;
fder = M4 * M1' .* (2*q^2 / sqrt(2*pi*v));
