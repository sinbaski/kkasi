clear all
close all

AR = [0:0.05:0.8, 0.9, 0.95:0.005:0.99];
j = 1;
tau = -log(2)/log(AR(j))
TailExponent = 2.9664;
% N = 50;
% q = 1600;
% T = N * q;

N = 4;
rho = 0;
q = 1500;
phi = 0;
% tau = -log(2)/log(0.99)

% hold on
% load(sprintf('GarchWishartT8e4Eig-%d.mat', j), 'ev');
% maxeig = max(ev);
% [x, y] = epdf(maxeig, 1, 0, 0.2, 400, 'b');
% ev = ev./mean(ev);
% param_ev = stblfit(ev, 'ecf')

% load(sprintf('GarchWishartT8e4Diag-%d.mat', j), 'Cii');
% Cii = reshape(Cii, prod(size(Cii)), 1);
% epdf(Cii, 1, 0, 0.1, 400, 'g');
% hold off
% Cii = Cii./mean(Cii);

%% plot the PDF of eigenvalues and the diagonal elements
% load(sprintf('GarchWishartT8e4Eig-%.3f.mat', AR(j)), 'ev');
% E = reshape(ev, prod(size(ev)), 1);
% [y_ev, x_ev] = ecdf(ev);
% [x_ev, y_ev] = epdf(E, 4, 0, 0.03, 500, 'b');
% params = stblfit(ev, 'ecf')
% y1 = stblcdf(x_ev, params(1), params(2), params(3), params(4));
% plot(log10(x_ev), log10(y_ev), 'b', ...
%      log10(x_ev(y1 > 0)), log10(y1(y1 > 0)), 'r');
% % plot((x_ev(y1 > 0)), (y1(y1 > 0)), 'r');

load(sprintf('GarchWishartN%dQ%dRho%.3fEig-%.3f.mat', ...
             N, q, rho, phi), 'ev');
% load(sprintf('GarchWishartT8e4Diag-%.3f.mat', AR(j)), 'Cii');
E = [...
    reshape(ev(1, :), 1, size(ev, 2));
    reshape(ev(2, :), 1, size(ev, 2));
    reshape(ev(3, :), 1, size(ev, 2));
    reshape(ev(4, :), 1, size(ev, 2));
    ];
F = reshape(E, 1, N*size(ev, 2));

% [y_cii, x_cii] = ecdf(F);
close all
hold on
[x_cii, y_cii] = epdf(F, 4, quantile(F, 0.2), quantile(F, 0.9), 400, 'g');
% param_Cii = stblfit(F, 'ecf');
y = stblpdf(x_cii, param_Cii(1), param_Cii(2), param_Cii(3), ...
            param_Cii(4));
plot(log10(x_cii), log10(y), 'r', 'LineWidth', 2);
grid on
hold off

% figure;
% load(sprintf('GarchWishartT8e4Offdiag-%.3f.mat', AR(j)), 'Cij');
% Cij = reshape(Cij, prod(size(Cij)), 1);
% % Cij = Cij(1:1e+6);
% [y_Cij, x_Cij] = ecdf(Cij);
% % [x_Cij, y_Cij] = epdf(Cij, 1, min(Cij), max(Cij), 500, 'g');
% param_Cij = stblfit(Cij, 'ecf')
% y1 = stblcdf(x_Cij, param_Cij(1), param_Cij(2), param_Cij(3), ...
%             param_Cij(4));
% hold on
% plot(x_Cij, log10(y_Cij), 'g', (x_Cij(y1 > 0)), log10(y1(y1 > 0)), 'm');
% hold off


% plot((x_cii(y>0)), (y(y > 0)), 'm');

% load(sprintf('GarchWishartT8e4Offdiag-%.3f.mat', AR(j)), 'Cij');
% Cij = reshape(Cij, prod(size(Cij)), 1);
% [y_cij, x_cij] = ecdf(Cij(1:1e+6));
% param_Cij = stblfit(Cij(1:1e+6), 'ecf')
% y1 = stblcdf(x_cij, param_Cij(1), param_Cij(2), param_Cij(3), param_Cij(4));

% figure;
% plot((x_cij), log10(y_cij), 'b', ...
%      (x_cij(y1 > 0)), log10(y1(y1 > 0)), 'g');






%% Fit the eigenvalue PDF to the custum PDF
% a = min(ev);
% [v I] = max(y_ev);
% l = x_ev(I);
% b = 1/(2/l - 1/a);
% b = 8.0e-3;
% v = (a + b + 2*sqrt(a*b))/4;
% q = (a + b -2*sqrt(a*b))/(a + b + 2*sqrt(a*b));
% plot(x_ev, MarcenkoPasturPDF(x_ev, [q, v]), 'g');

% thepdf = @(x, q, v) ...
%          [MarcenkoPasturPDF(x(x < b), [q, v]);...
%           stblpdf(x(x >= b), params(1), params(2), ...
%                   params(3), params(4))];
% [mp, bounds] = mle(ev, 'pdf', thepdf, 'start', [q, v]);

% params = mle(ev, 'pdf', @(x, u) exp(-u) * stblpdf(x, param_Cii(1),...
%                                                   1, ...
%                                                   param_Cii(3),...
%                                                   param_Cii(4)), ...
%              'start', 1, 'lowerbound', 0);
% plot(x_ev, thepdf(x_ev', mp(1), mp(2)), 'g');
% hold off
% legend('Eigenvalues', 'Cii', ...
%        sprintf('S(%.2f, %.2f, %.2e, %.2e)', param_Cii));
% title(sprintf('Eigenvalues and Diagonal Elements Distribution. \\tau = %.2f', ...
%               tau), 'Fontsize', 10);
% grid on

%% Plot the non-diagonal elements' PDF
% figure;
% load(sprintf('GarchWishartT8e4Offdiag-%.3f.mat', AR(j)), 'Cij');
% Cij = reshape(Cij, prod(size(Cij)), 1);
% param_Cij = stblfit(Cij, 'ecf')
% hold on
% % x = epdf(Cij, 1, -3e-4, 3e-4, 1000, 'b');
% x = epdf(Cij, 1, -0.01, 0.01, 1000, 'g');
% y = stblpdf(x, param_Cij(1), param_Cij(2), param_Cij(3), ...
%             param_Cij(4));
% plot(x, y, 'r');
% hold off
% title(sprintf('Cij distributions. \\tau=%.2f', tau));
% legend('Cij', sprintf('S(%.2f, %.2f, %.2e, %.2e)', param_Cij));
% grid on



% x = 0:1e-3:0.1;
% y1 = LevyPDF(x, 1, 1, 1, 1);
% y2 = stblpdf(x, 1, 1, 1, 1);
% plot(x, y1, x, y2);

% param = mle(X, 'pdf', @(x, beta, gam) stblpdf(x, alfa/2, beta, gam, 0),...
%             'start', [1, (2*pi)^(1/alfa) * (c0_stem * c1)^(2/alfa)]);

% x = min(Cii) : 2e-4:0.1;
% y = hist(Cii, x) ./ length(Cii) ./ (x(2) - x(1));
% y1 = stblpdf(x, param(1), param(2), param(3), param(4));
% plot(x(1:end-1), y(1:end-1), x, y1);



% a0 = 2.3e-6;
% a1 = 0.90;
% b1 = 0.09;
% a1 = 0.1;
% b1 = 0.8;
% DoF = 3;

% T = 1e+5;
% dist = struct('name', 'Garch1_1', 'prmt', [a0, a1, b1], 'distr', ...
%               struct('Name', 't', 'DoF', 3), 'TailExponent', 0.7744);
% % dist = struct('name', 'Garch1_1', 'prmt', [a0, a1, b1], ...
% %               'distr', struct('Name', 'Gaussian'));
% R = gen_ret_mtx(1, T, dist, 0);
% [P, x] = ecdf(R);
% P = 1 - P;
% I = x > exp(-4.5) & x < exp(-2);
% plot(log(x(I)), log(P(I)));
% grid on
% param = polyfit(log(x(I)), log(P(I)), 1)

%% Calculate the tail exponent of the returns
% T = 1e+6;
% Z = randn(1, T);
% % Z = trnd(3, 1, T);
% % Z = Z ./ sqrt(3);
% A = a1.* Z.^2 + b1;
% s = mean(log(A))
% if s < 0
%     f = @(k) mean(A.^(k/2)) - 1;
%     fsolve(f, 3)
% end

% I = find(y > 4*s & y < 10*s & P > 0);
% plot(log(y(I)), log(P(I)), 'o-');
% grid on
% coef = polyfit(log(y(I)), log(P(I)), 1)

%%
% log P = log(c) - 2 log|x|
%
% c = exp(-9.28) = 9.33e-5;
% Cn = 2*T*c / pi;


