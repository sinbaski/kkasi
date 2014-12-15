clear all
close all

VaR = 219.27;
S = [
    1, 1/4, 1/3
    1/4, 1, -1/5
    1/3, -1/5, 1]*1.0e-3;
Mu = 1.0e-3 * [1/2, 1, 1]';
% a = 1.0e+3;
% b = (-VaR/a - sum(Mu))/()


Nu = [0; 0; -219.27e-3];
num = 5.0e+3;
R = mvnrnd(Nu, S, num);
W = ones(1, 3);
Ln1 = sort(-sum(diag([1000, 1000, 1000]) * R'));

hold on; plot(R(1, :), -VaR*1.0e-3 - R(1, :), 'r-');







