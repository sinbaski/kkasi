close all
% f = @(p, x) p(1) .*  tan(pi*x/2) + p(2);
f = @(p, x) p(1) * sin(pi/2./(1 - x))./ (1 - x);
%% alpha in the magnitude
% p = [ph_params(2, 11) - ph_params(2, 1), ph_params(2, 1)];
%% alpha in the phase
% p = [ph_params(2, 11) - ph_params(2, 1), ph_params(2, 1)];
% p = [ph_params(2, 1) - ph_params(2, 11), ph_params(2, 11)];
% p = [ph_params(2, 1)];
% P = lsqcurvefit(f, p, AR(collection), ph_params(2, collection));
% x = linspace(0, 0.99, 500);
% figure;
% plot(AR(collection), ph_params(3, collection), '-o');
x = AR(collection);
% y = ph_params(3, collection);
% plot(pi / 2 ./ (1 - x), y .* (1 - x));
plot(x, ph_params(3, collection), 'o-');
hold on
% plot(x, mag_params(1, collection), 'ro-');
grid on
hold off
