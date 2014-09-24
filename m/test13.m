clear all
close all

X = [1:0.001:2]' * 1.0e-3;
Y = [-0.5:0.01:0.5]'.*1.0e-4;

v = 1;

M = NaN(length(X), length(Y));
for n = 1 : length(Y)
    M(:, n) = LognormalSigmaM(X + i*Y(n), v);
end

subplot(1, 2, 1);
contour(X, Y, real(M'), 10);
ylim([-2, 2]*1.0e-6);
xlabel('Re Z');
ylabel('Im Z');
colorbar('location', 'EastOutside');
title('Re M');

subplot(1, 2, 2);
contour(X, Y, imag(M'), 10);
ylim([-0.5, 0.5]*1.0e-5);
xlabel('Re Z');
ylabel('Im Z');
colorbar('location', 'EastOutside');
title('Im M');
