clear all
close all
n = 500;
X = linspace(-8, 1, n)';
Y = linspace(-4, -1.0e-3, n)';

Bre = NaN(n, n);
Bim = NaN(n, n);
for k = 1 : n
    M = LognormalBlue(X(k) + i*Y, 0.5, 0.5);
    Bre(:, k) = M(:, 1);
    Bim(:, k) = M(:, 2);
end
C = contour(X, Y, Bim, [0, 0]);

subplot(1, 2, 1);
contourf(X, Y, Bre, [-1:0.1:1]);
colorbar('location', 'EastOutside');
hold on
plot(C(1, 2:end), C(2, 2:end), 'k-', 'Linewidth', 2);
title('Re B');
xlabel('Re G');
ylabel('Im G');

subplot(1, 2, 2);
contourf(X, Y, Bim, [-1:0.1:1]);
colorbar('location', 'EastOutside');
title('Im B');
xlabel('Re G');
ylabel('Im G');




% A = find(abs(Bim) < 1.0e-4);
% C = NaN(length(A), 2);
% for k = 1 : length(A)
%     C(k, 1) = ceil(A(k)/n);
%     C(k, 2) = mod(A(k), n);
%     if (C(k, 2) == 0)
%         C(k, 2) = n;
%     end
% end
% hold on
% plot(X(C(:, 1)), Y(C(:, 2)), 'k.', 'Linewidth', 2);
