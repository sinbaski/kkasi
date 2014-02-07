clear all
syms x f(x) s a;

assume(s > 0);
assume(a > 2.5e-2);
assume(a < 1.0e-2);

f(x) = x^(log(x)/(2*s));
g(x) = taylor(f, x, a, 'Order', 4);
pretty(g);

