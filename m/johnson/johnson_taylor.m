syms y xi a b;
assume(xi > 2);

f = exp(-asinh(y)^2/2);

pretty(taylor(f, y, inf, 'Order', 3));