clear all
close all

S = [
    1, 1/4, 1/3
    1/4, 1, -1/5
    1/3, -1/5, 1]*1.0e-3;

VaR = 219.27;
Nu = [0; 0; -VaR*1.0e-3];
Mu = 1.0e-3 * [1/2, 1, 1]';
xi = S^(-1) * (Nu - Mu);
C = exp(-(2*xi'*Mu + xi'*S^2*xi)./2);
D = 1/C;
num = 1.0e+4;
R = mvnrnd(Nu, S, num)';

U = linspace(VaR-10, VaR+10, 1.0e+2)';
% U = [VaR-10, VaR, VaR+10];
Fbar = NaN(length(U), 1);
L = -1.0e+3*sum(R);
for k = 1 : length(U)
    I = find(L > U(k));
    Fbar(k) = sum(D*exp(-xi' * R(:, I)))/num;
end
plot(U, Fbar);
