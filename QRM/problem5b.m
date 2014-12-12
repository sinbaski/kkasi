clear all
close all

S = [
    1, 1/4, 1/3
    1/4, 1, -1/5
    1/3, -1/5, 1]*1.0e-3;

Mu = [0; 0; -219.27e-3];
num = 1.0e+5;
R = mvnrnd(Mu, S, num);
W = ones(1, 3);
Ln1 = sort(-W * diag([1000, 1000, 1000]) * R');
VaR = 219.27;







