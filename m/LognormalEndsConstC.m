function C = LognormalEndsConstC(v, q, a, n)
L = 4000;
s = linspace(-6*sqrt(v), 6*sqrt(v), L);
M1 = exp(-s.^2./2./v)*(s(2) - s(1))/sqrt(2*pi*v);
M2 = 1./(exp(-2.*s) - q.*a).^n;
C = M2 * M1';

