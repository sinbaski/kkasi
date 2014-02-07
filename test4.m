clear all
syms x y z p f(x,y);

assume(p < 1);
assumeAlso(p > -1);
assume(z, 'real');


% f(x,y) = 1/(2*pi)/sqrt(1 - p^2) * exp(-(x^2 + 2*p*x*y + y^2)/2);
f(x) = exp(-)

g = int(f(x,y), x, 0, z/y);
