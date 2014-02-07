a = -25.59;
b = 46.54;
A = 2.0e-4;
B = -2.3e-5;
C = 2.6e-4;
D = -3.261e-3;
E = -4.667e-2;
F = 69.22;

f = @(t) D*t.^2 + E*t + F;
q = @(t) A*t.^2 + B*t + C;
phi = @(t) 2.^(-1./t);
k = @(t) a*phi(t) + b;

t = linspace(min(tau), max(tau), 200);

M1 = @(t) k(t).*q(t) - f(t);
M2 = @(t) k(t).*q(t).^2;
M3 = @(t) 2./sqrt(k(t));


subplot(1, 2, 1);
plot(tau, tw_moments(:, 1), 'bx', t, M1(t), 'r');

subplot(1, 2, 2);
plot(tau, tw_moments(:, 2), 'bx', t, M2(t), 'r');

% subplot(1, 3, 3);
% plot(tau, tw_moments(:, 3), 'bx', t, M3(t), 'r');
