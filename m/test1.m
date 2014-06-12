clear all
dist = struct('name', 'Gaussian');
k = NaN(18*19, 1);
c = 1;
a0 = 0.01;
for a1 = 0.05:0.05:0.9
    n = 10;
    db = (1 - a1)/n;
    for b1 = 0:0.05:0.9
        param = [0.01, a1, b1];
        k(c) = garch_tail_exponent(dist, param, 1e+6);
        c = c + 1;
    end
end
save('garch_tail_exponent.txt', 'k', '-ascii');


