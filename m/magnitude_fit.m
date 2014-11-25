
y1 = params(:, 2)./params(:, 1);
x1 = AR(collection);
tau1 = -log(2) ./ log(x1);

% figure;
% plot(x1(1:17), y1(1:17), 'o');

a = 0.3180;
b = -7.165;
% hold on
plot(x1(1:17), atan((y1(1:17) - b)./a), 'o');

f = @(p, x) p .* x;
P = lsqcurvefit(f, pi/2, x1(1:17), atan((y1(1:17) - b)./a)');
hold on
plot(x1(1:17), f(P, x1(1:17)), 'g');

