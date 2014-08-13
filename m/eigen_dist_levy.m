clear all
close all
% figure;
T = 8e+4;
N = 50;
fold = 4000;

%% eigen values generated from the volvo stochastic log-volatility
%% model.

% spec = cellstr(['b  '; 'g  '; 'r  '; 'k  ']);
spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
                'b: '; 'c: '; 'g: '; 'm: '; 'r: '; 'k: ';...
                'b-.'; 'c-.'; 'g-.'; 'm-.'; 'r-.'; 'k-.';...
                ]);

%% Approximation of the Tracy-Widom with a Gamma dist.
% a1 = -1/2;
% a2 = a1;
% n = N;
% p = T;
% mu = (sqrt(n+a1) + sqrt(p+a2))^2;
% sigma = sqrt(mu) * (1/sqrt(n+a1) + 1/sqrt(p+a2))^(1/3);

% hold on;
% sample_moments = NaN(length(AR), 3);
% dev_moments = NaN(length(AR), 3);
% gam_moments = NaN(length(AR), 3);
% tw_moments = NaN(length(AR), 3);
% alpha = NaN(length(AR), 1);

c = 1;
% maxeig = NaN(length(collection), fold);
% mineig = NaN(length(collection), fold);
% mag_params = NaN(2, length(collection));
% ph_params = NaN(3, length(collection));
% Cii_params = NaN(4, 4);
% ev_params = NaN(9, 4);
% merge_point = NaN(9, 4);
N = 50;
q = 24;

% alpha = 2.17
% garch_params = [
%     0.1068,  0.8923
%     % 0.2009, 0.7962
%     0.2991,  0.6946
%     % 0.4001, 0.5889
%     0.4996,  0.4835
%     % 0.5995, 0.3786
%     0.6992,  0.2717
%     % 0.8009, 0.1608
%     % 0.8970,  0.0541
%         ];

a = (1 - sqrt(1/q))^2;
b = (1 + sqrt(1/q))^2;

% alpha = 3.00
% garch_params = [
%     0.5000, 0.3874 % a = 2.9936, tau = 0.9999
%     0.6031, 0.2259 % a = 3.0088, tau = 1.9996    
%     0.7156, 0.0289 % a = 3.0099, tau = 2.2151
%                ];

% alpha = 4.00
% garch_params = [
%     0.1012, 0.8885 % , a = 4.0097
%     0.2996, 0.6060 % , a = 4.0027, tau = -4.1641    
%     0.5004, 0.2082 % , a = 3.9908, tau = 1.0023
%          ];

% alpha = 4.00
% garch_params = [
%     0.5004, 0.2082 % a = 3.9908, tau = 1.0023
%     0.5491, 0.0784 % a = 4.0086, tau = 1.2009    
%                ];

garch_params = [
    0.4996, 0.4835, 2.17 % , tau = 0.9530
    0.5000, 0.3874, 2.9936 % , tau = 0.9999
    0.5004, 0.2082, 3.9908 %, tau = 1.0023    
               ];

a1 = -1/2;
a2 = a1;
n = N;
p = T;
mu = (sqrt(n+a1) + sqrt(p+a2))^2;
sigma = sqrt(mu) * (1/sqrt(n+a1) + 1/sqrt(p+a2))^(1/3);

% gam_moments = NaN(length(garch_params), 3);
% tw_moments = NaN(length(garch_params), 3);
% alpha = NaN(length(garch_params), 1);
% k = NaN(length(garch_params), 1);
% theta = NaN(length(garch_params), 1);

% var_lambda = NaN(length(garch_params), 1);

texts = {};
hold on
while c <= size(garch_params, 1)
    % load(sprintf('../data/GarchWishartOrdN%dQ%dA%.4fB%.4fEig.mat',...
    %              N, q, garch_params(c, 1), garch_params(c, 2)), 'ev');
    % maxeig = max(ev);
    % maxeig1 = maxeig;
    % % % maxeig1 = (maxeig - mean(maxeig))./std(maxeig);
    % [x1, y1] = epdf(maxeig1, 1, min(maxeig1), max(maxeig1), 100, ...
    %                 '');
    % stairs(x1, y1, spec{c}, 'LineWidth', 2);
    % fprintf('<lam1> = %.4f, var(lam1) = %.4f\n', mean(maxeig1), var(maxeig1));
    % % % [y1, x1] = ecdf(maxeig1);
    % % % plot(log10(x1(x1 > 1)), log10(1-y1(x1 > 1)), spec{c});
    % texts{c} = sprintf('\\alpha = %.2f', exponent(c));
    % c = c + 1;
    
    % maxeig1 = (maxeig - mu)/sigma;
    % tw_moments(c, :) = [mean(maxeig1), var(maxeig1), ...
    %                     skewness(maxeig1)];
    % moments = tw_moments(c, :);

    % k(c) = 4/(moments(3)^2);
    % theta(c) = sqrt(moments(2)) * moments(3)/2;
    % alpha(c) = k(c)*theta(c) - moments(1);
    % maxeig1 = maxeig1 + alpha(c);
    
    % [y, x] = ecdf(maxeig1);
    % y1 = cdf('Gamma', x, k(c), theta(c));
    % var_lambda(c) = max(abs(y - y1));
    % %    plot((x), y, spec{c}, x, y1, 'k', 'LineWidth', 2);
    % c = c + 1;
    
    % [y1, x1] = ecdf(Cii);
    % y2 = stblcdf(x1, Cii_params(1), Cii_params(2), Cii_params(3), Cii_params(4));
    % [y3, x3] = ecdf(ev);
    % plot(log10(x3), log10(y3), log10(x1), log10(y1), log10(x1), ...
    %      log10(y2));
    % legend('Eigenvalue', 'Cii', sprintf('S(%.2f, %.2f, %.2e, %.2e)', ...
    %                                     Cii_params(1), Cii_params(2), ...
    %                                     Cii_params(3), Cii_params(4))); 
    % xlim([-2.3, -1.9]);
    % ylim([-3.5, 0]);
    % grid on
    % title(sprintf('N=250, T=%d', N*q));
    % maxeig = max(ev);
    % [y, x] = ecdf(maxeig);
    % plot(log(x), log(y), spec{c});
    % c = c + 1;

    %% Plot the PDF of the maximum eigenvalue / diagonal elements
    % load(sprintf('GarchWishartWrongT8e4Eig-%.3f.mat', AR(j)), 'ev');
    % ev = reshape(ev, prod(size(ev)), 1);
    % [x, y] = epdf(ev, 1, 0, 0.05, 400, spec{c});
    % params = stblfit(ev, 'ecf');
    % y1 = stblpdf(x, params(1), params(2), params(3), params(4));
    % plot(x, y, x, y1);
    % c = c + 1;

    %% spectral distribution
    a0 = 1 - sum(garch_params(c, :));
    a1 = garch_params(c, 1);
    b1 = garch_params(c, 2);
    fprintf('a1 = %.2f, b1 = %.2f\n', a1, b1);
    
    load(sprintf('../data/GarchWishartOrdN%dQ%dA%.4fB%.4fEig.mat',...
                 N, q, a1, b1), 'ev');
    ev = reshape(ev, prod(size(ev)), 1);
    fprintf('<eigenvalue> = %.4f, var(eigenvalue) = %.4f\n', mean(ev), var(ev));
    [x1, y1] = epdf(ev, 4, min(ev), max(ev), 400, '');
    % metric(c) = HellingerDistance(...
    %     @(x) y1, ...
    %     @(x) MarcenkoPasturPDF(x, [1/q, 1]), x1);

    plot(log10(x1), log10(y1), spec{c}, 'LineWidth', 2);
    % [y1, x1] = ecdf(ev);
    % plot((x1(x1 > 1.5)), log10(y1(x1 > 1.5)), spec{c}, ...
    %      'LineWidth', 2);
    grid on
    texts{c} = sprintf('\\rho_n = %.4f^{n-1} %.2f, \\alpha=%.2f',...
                      sum(garch_params(c, 1:2)), garch_params(c,...
                      1), garch_params(c, 3));
    % xlabel('Eigenvalue', 'Fontsize', 14);
    % ylabel('PDF of Eigenvalue', 'Fontsize', 14);
    c = c + 1;
    
    % m = N*(N-1)/2;
    % % m = N;
    % Cij = NaN(1, m * fold);
    % load(sprintf('../data/GarchWishartOrdN%dQ%dA%.4fB%.4fC.mat', N, ...
    %              q, a1, b1), 'C');
    % for n = 0:size(C, 3)-1
    %     D = C(:, :, n+1);
    %     Cij(n*m+1 : (n+1)*m) = D(logical(triu(ones(N), 1)));
    %     % Cij(n*m+1 : (n+1)*m) = D(logical(eye(N)));
    % end
    % % [x1, y1] = epdf(Cij, 1, -0.15, 0.15, 1000, '');
    % % plot(x1, (y1), spec{c}, 'LineWidth', 2);
    % [y1, x1] = ecdf(Cij);
    % plot(log10(x1(x1 > 0)), log10(1 - y1(x1 > 0)), spec{c}, ...
    %      'LineWidth', 2);
    % grid on
    % texts{c} = sprintf('\\rho_n = %.4f^{n-1} %.2f',...
    %                   sum(garch_params(c, :)), garch_params(c, 1));
    % xlabel('Cij', 'Fontsize', 14);
    % ylabel('CDF of Cij', 'Fontsize', 14);
    % c = c + 1;
    
    % p2 = stblfit(Cii, 'ecf');
    % y3 = stblcdf(x2, p2(1), p2(2), p2(3), p2(4));
    % Cii_params(:, c) = p2;

    % plot(log10(x1), log10(y1), 'b', log10(x2), log10(y2), 'g', ...
    %      log10(x2), log10(y3), 'r');
    % title(sprintf('N=250, T=%d', T));

    % y1 = stblpdf(x1, p1(1), p1(2), p1(3), p1(4));
    % plot(log10(x1(y1 > 0)), log10(y1(y1 > 0)), 'r');
    % y2 = stblpdf(x2, p2(1), p2(2), p2(3), p2(4));

    % plot(log10(x2(y2 > 0)), log10(y2(y2 > 0)), 'm');
    % grid on
    % hold off
    % c = c + 1;
    
    % figure;
    % plot(log(x), log(y), 'b', log(x), log(y1), 'r');
    % legend(sprintf('\\phi=%f', (j-1)*0.05),...
    %        sprintf('\\Gamma(%f, %f)', k, theta), 'Location', 'Northeast');
    % xlabel('(\lambda_1 - \mu)/\sigma');
    % ylabel('PDF of (\lambda_1 - \mu)/\sigma');
    
    %% plot the CDF of the original
    % [y, x] = ecdf(maxeig);
    % plot(x, y, spec{j});

    %% plot the CDF of the transformed
    % [y, x] = ecdf(maxeig1);
    % y1 = cdf('Gamma', x, k, theta);
    % subplot(2, 2, (j-1)/3+1);
    % plot(log(x), log(y), 'b', log(x), log(y1), 'r');
    % legend(sprintf('\\phi=%f', (j-1)*0.05),...
    %        sprintf('\\Gamma(%f, %f)', k, theta), 'Location', 'Southeast');
    % xlabel('ln((\lambda_1 - \mu_{NT})/\sigma_{NT})');
    % ylabel('ln(F(\lambda_1))');
    
    %% plot CDF of maxeig - (1 + sqrt(q))^2
    % dev = maxeig - (1 + sqrt(q))^2 / (1 - AR(j)^2);
    % dev_moments(j, :) = [mean(dev), std(dev), skewness(dev)];
    % [y, x] = ecdf(dev);
    % y1 = cdf('Normal', x, mean(dev), std(dev));
    % plot(x, y, x, y1, spec{j});
    
    
    %% Distance between Cii and ev
    % subplot(3, 3, c);
    % load(sprintf('GarchWishartN%dQ%dEig-%.3f.mat', N, q, 0), 'ev');
    % ev = sort(reshape(ev, prod(size(ev)), 1));
    
    % load(sprintf('GarchWishartN%dQ%dDiag-%.3f.mat', N, q, 0), 'Cii');
    % Cii = sort(reshape(Cii, prod(size(Cii)), 1));

    % distance = ev - Cii;
    % I = abs(distance) ./ min([Cii'; ev'])' < 1.0e-3;
    % I(1:16) = 0;
    % J = I;
    % for n = 1:16
    %     J = circshift(J, -1);
    %     I = I & circshift(J, 1);
    % end
    % a = min(ev(I));
    
    % hold on
    % epdf(ev, 4, 0, 2.0e-2, 1000, 'b');
    % epdf(Cii, 4, 0, 2.0e-2, 1000, 'g');

    % epdf(ev, 1, a, 2.0e-2, 500, 'r');
    % epdf(Cii, 1, a, 2.0e-2, 500, 'm');
    % hold off
    % grid on
    % merge_point(c, (q-12)/4+1) = a - mean(ev);
    % c = c + 1;
    
    % fprintf('Done with #%d.\n', c);
end

% c = 1;
% while c <= size(garch_params, 1)
%     texts{c} = sprintf('[%.4f, %.4f]', garch_params(c, :));
%     c = c + 1;
% end
% x2 = linspace(a, b, 10000);
% y2 = MarcenkoPasturPDF(x2, [1/q, 1]);
% plot((x2), log10(y2), 'k', 'LineWidth', 2);
hold off
% plot(garch_params(:, 1), (var_lambda), 'b*-');
grid on

set(gca, 'Fontsize', 14);
legend(texts);
%  legend('\tau_0 = 1.00', '\tau_0 = 2.00', '\tau_0 = 2.22');
% title('P(|r| > x) ~ 1/x^{2.17}', 'Fontsize', 14);
xlabel('log(eigenvalues)', 'Fontsize', 14);
ylabel('log(spectral density)', 'Fontsize', 14);

% if exist('garch_eig_fourier_coef_T8e4.mat', 'file') ~= 2
%     save('garch_eig_fourier_coef_T8e4q.mat', 'C');
% end

%% Plot p against phi
% figure;
% subplot(1, 3, 1);
% plot(AR(collection), p(1, collection));

% subplot(1, 3, 2);
% plot(AR(collection), p(2, collection));

% subplot(1, 3, 3);
% plot(AR(collection), abs(p(3, collection)));

% grid on

% m1 = gam_moments(:, 1);
% m2 = gam_moments(:, 2);
% m3 = gam_moments(:, 3);

% k = gam_moments(:, 1).^2 ./ gam_moments(:, 2);
% theta = gam_moments(:, 2) ./ gam_moments(:, 1);

% p1 = polyfit(tau, k, 1);
% p2 = polyfit(tau, theta, 1);
% plot(tau, k, 'bx', tau, polyval(p1, tau), 'r');
% hold on
% plot(tau, theta, 'co', tau, polyval(p2, tau), 'm');
% hold off
% grid on
% legend('k', sprintf('k=%.1f \\tau %+.1f', p1(1), p1(2)), '\theta', ...
%        sprintf('\\theta=%.2f \\tau %+.1f', p2(1), p2(2)));
% xlabel('\tau');


% plot(tau, log(m1-min(m1)), tau, log(m2), tau, log(m3));
% plot(log(AR), log(m1), 'bx-');
% P = polyfit(tau, m1', 2);
% plot(tau, m1, 'bo-', tau, polyval(P, tau), 'r');
% plot(m1, m2, 'bx-');
% xlabel('\tau');
% ylabel('ln(<\lambda_{max}>)')
% load('GaussianWishartEig-1.mat', 'ev');
% ev = reshape(ev, 1, N*fold);
% lag = 1;
% x = linspace(min(ev)*0.8, max(ev), 100);
% y = hist(ev, x) / length(ev) / (x(2) - x(1));
% x = tsmovavg(x, 's', lag);
% x = x(lag:end);
% y = tsmovavg(y, 's', lag);
% y = y(lag:end);
% y1 = MarcenkoPasturPDF([N/T, 1], x);
% figure;
% plot(x(51:end), y(51:end), x(51:end), y1(51:end));

% legend('$\phi=0$', '$\phi=0.5$', '$\phi=0.8$', '$\phi=0.955$', '$\phi=0.97$', ...
%        'Location', 'Northeast');
% h = legend;
% set(h, 'interpreter', 'latex', 'fontsize', 14);

