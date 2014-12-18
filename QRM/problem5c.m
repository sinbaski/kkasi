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
R = mvnrnd(Nu, S, num)';
L = -1.0e+3*sum(R);

Q = linspace(VaR-5, VaR+5, 4.0e+2)';
Fbar = NaN(length(Q), 1);
for k = 1 : length(Q)
    I = find(L > Q(k));
    Fbar(k) = sum(D*exp(-xi' * R(:, I)))/num;
end
plot(Q, Fbar, '.');
xlabel('x', 'Fontsize', 14);
ylabel('$\bar{F}(x)$', 'Interpreter', 'latex', 'Fontsize', 14);
grid on

k = min(find(Fbar < 1.0e-4));
fprintf('Estimated VaR is %.4f\n',Q(k));
