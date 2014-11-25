function density = HestonStub(param, x)
model = heston(param(1), param(2), param(3), param(4), 'Correlation', ...
               [1, 0; param(5), sqrt(1 - param(5)^2)]);
t = 30;
density = HestonRetDistr(model, t, x);
