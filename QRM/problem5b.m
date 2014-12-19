clear all
close all

S = [
    1, 1/4, 1/3
    1/4, 1, -1/5
    1/3, -1/5, 1]*1.0e-3;

% VaR = 219.27;
VaR = 162.5;
% Nu = -[1; 1; 1] * VaR*1.0e-3/3;
% Nu = -[0; 0; 1] * VaR*1.0e-3;
Mu = 1.0e-3 * [1/2, 1, 1]';
% xi = S^(-1) * (Nu - Mu);
xi = -(VaR*1.0e-3 + sum(Mu))/sum(sum(S)) * ones(3, 1);
Nu = Mu + S*xi;
C = exp(-(2*xi'*Mu + xi'*S^2*xi)./2);
D = 1/C;
num = 1.4e+4;
rounds = 2000;

prob = NaN(1, rounds);
for k = 1 : rounds
    R = mvnrnd(Nu, S, num)';
    L = -1.0e+3*sum(R);
    I = find(L > VaR);
    prob(k) = sum(D*exp(-xi' * R(:, I)))/num;
end
prob = sort(prob);

m = rounds;
alfa = 1 - 1.0e-4;
accum = 0;
bit = 1.0e-2;
j = 0;
while accum < bit/2
    A = [m : -1 : m-j+1] .* alfa;
    B = (1 - alfa) ./ [1:m-j];
    accum = accum + prod(A)*prod(B);
    j = j + 1;
end
j = j - 1;
l = j + 1;
k = m - j;

fprintf('P(L > VaR) in [%.2e, %.2e] with %.2f confidence\n', prob(l), ...
        prob(k), 1- bit);






