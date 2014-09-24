clear all
close all
v = 1.0e-4;
q = 0.1;

a = (1 - sqrt(q))^2;
b = (1 + sqrt(q))^2;

% [X, Y] = meshgrid([-0.35:0.005:0.35], [-0.5:0.01:-0.1]);
% Z = X + i.*Y;
% F = NaN(size(Z, 1), size(Z, 2), 2);
% for n = 1:size(Z, 2)
%     B = LognormalBlue(Z(:, n), v, q);
%     F(:, n, 1) = B(:, 1);
%     F(:, n, 2) = B(:, 2);
% end
% subplot(1, 2, 1);
% contourf(X, Y, F(:, :, 1), 20);
% xlabel('Re G');
% ylabel('Im G');
% colorbar('location', 'EastOutside');
% title('Re B');

% subplot(1, 2, 2);
% contourf(X, Y, F(:, :, 2), 20);
% xlabel('Re G');
% ylabel('Im G');
% colorbar('location', 'EastOutside');
% title('Im B');

% [min(min(F(:, :, 2))), max(max(F(:, :, 2)))]

% g0 = LognormalGreen((a+b)/2, v, q);
z0 = 0.5;
flag = 1;
while flag
    g0 = LognormalGreen(z0, v, q);
    if (abs(imag(g0)) > 1.0e-3)
        flag = 0;
    else
        z0 = z0*1.05;
    end
end
[end1, end2] = LognormalEnds(v, q, g0);
Z1 = linspace(end1, end2, 500)';
G1 = LognormalGreen(Z1, v, q);

Z2 = linspace(a, b, 500)';
P2 = MarcenkoPasturPDF(Z2, [q, 1])';

plot((Z1), log10(-imag(G1)./pi), 'b', (Z2), log10(P2), ...
     'r', 'LineWidth', 2);
grid on
legend('Lognormal SV', 'MP');
ylabel('Spectral density function');

% Z3 = linspace(end1, b, 500)';
% G3 = LognormalGreen(Z3, v, q);

% plot((Z3), log10(-imag(G3)./pi), 'b', (Z2), log10(P2), ...
%      'r', 'LineWidth', 2);
% xlim([min(a, end1), b]);
% grid on
