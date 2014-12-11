clear all
close all

Mu = [1/2, 1, 1] * 1.0e-3;

S = [
    1, 1/4, 1/3
    1/4, 1, -1/5
    1/3, -1/5, 1];

alfa = 1- 1.0e-4;
bita = 1 - 0.99;

num = 1.0e+5;
R = mvnrnd(Mu, S, num);
W = ones(1, 3);
Vn1 = W * diag([1000, 1000, 1000]) * exp(R');
Ln1 = sort(3000 - Vn1, 'descend');
VaR = Ln1(ceil(num*(1-alfa)));
[lb, ub] = VaRBounds(1 - 1.0e-4, 0.99, num);






