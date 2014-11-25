clear all
syms a s p v m x

assume(p < 1);
assumeAlso(p > -1);

assume(s > 0);

f = int(exp(-a^2/2) * (a * p  * exp(a * s + v) + m)^4, a, -inf, inf);
pretty(f);
