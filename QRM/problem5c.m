clear all
close all

S = [
    1, 1/4, 1/3
    1/4, 1, -1/5
    1/3, -1/5, 1]*1.0e-3;

% VaR = 219.27;
VaR = 167.5;
Nu = -[1; 1; 1] * VaR*1.0e-3/3;
% Nu = -[0; 0; 1] * VaR*1.0e-3;
Mu = 1.0e-3 * [1/2, 1, 1]';
xi = S^(-1) * (Nu - Mu);
C = exp(-(2*xi'*Mu + xi'*S^2*xi)./2);
D = 1/C;
num = 1.4e+5;
R = mvnrnd(Nu, S, num)';
L = -1.0e+3*sum(R);

U = linspace(VaR-20, VaR+20, 4.0e+2)';
Fbar = NaN(length(U), 1);
for k = 1 : length(U)
    I = find(L > U(k));
    Fbar(k) = sum(D*exp(-xi' * R(:, I)))/num;
end
plot(U, Fbar, '.');
grid on
