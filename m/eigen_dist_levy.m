clear all
close all

T = 8e+4;
N = 50;
fold = 2e3;

%% eigen values generated from the volvo stochastic log-volatility
%% model.

spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
                'b: '; 'c: '; 'g: '; 'm: '; 'r: '; 'k: ';...
                'b-.'; 'c-.'; 'g-.'; 'm-.'; 'r-.'; 'k-.';...
                ]);
% for i = [1:5]
%     load(sprintf('./eigen-%d.mat', i));
%     ev = reshape(ev, N, fold);
%     maxeig = max(ev);
%     [y, x] = ecdf(maxeig);
%     y1 = diff(y) ./ diff(x);
%     x1 = (x(1:end-1) + x(2:end))/2;
%     plot(log(x1), log(y1), spec(i));
%     hold on
% end
% hold off
% grid on


%% White Wishart matrix eigen values
% AR = [0:0.05:0.8, 0.9, 0.95:0.005:0.99];
AR = 0:0.1:0.9;
tau = -log(2)./log(AR);
% dist = struct('name', 'Cauchy', 'prmt', 1);
dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.15, 0.84], 'distr', ...
              'Gaussian', 'TailExponent', 2.9664);
% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.1, 0.80], 'distr', ...
%               struct('Name', 't', 'DoF', 3), 'TailExponent', 0.7744);
% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.1, 0.725], 'distr', ...
%               struct('Name', 't', 'DoF', 3), 'TailExponent', 1.8187);
% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.15, 0.84], 'distr', ...
%               struct('Name', 't', 'DoF', 60), 'TailExponent', 2.4477);
% for j = AR
%     if (exist(sprintf('GarchWishartT8e4Eig-%.3f.mat', j), 'file') == 2)
%         continue;
%     end
%     ev = NaN(N, fold);
%     Cij = NaN(N*(N-1)/2, fold);
%     Cii = NaN(N, fold);
%     for i = 1:fold
%         R = gen_ret_mtx(N, T, dist, j);
%         if strcmp(dist.name, 'Cauchy') == 1
%             R = R ./ dist.prmt ./ T;
%         elseif strcmp(dist.name, 'Garch1_1') == 1
%             % R = diag(1./std(R')) * R;
%             R = R ./ T^(1/dist.TailExponent);
%         end
%         C = R*R';
%         Cij(:, i) = C(logical(triu(ones(N), 1)));
%         Cii(:, i) = C(logical(eye(N)));
%         ev(:, i) = eig(C);
%     end
%     save(sprintf('GarchWishartT8e4Eig-%.3f.mat', j), 'ev');
%     save(sprintf('GarchWishartT8e4Offdiag-%.3f.mat', j), 'Cij');
%     save(sprintf('GarchWishartT8e4Diag-%.3f.mat', j), 'Cii');
% end

% [y, x] = ecdf(maxeig);
% p = RealWishartCDF(x(1:10:end), N, T);

% plot(x, y, x(1:10:end), p);

%% Approximation of the Tracy-Widom with a Gamma dist.
% a1 = -1/2;
% a2 = a1;
% n = N;
% p = T;
% mu = (sqrt(n+a1) + sqrt(p+a2))^2;
% sigma = sqrt(mu) * (1/sqrt(n+a1) + 1/sqrt(p+a2))^(1/3);

% hold on;
% avgs = NaN(1, length(AR));

% sample_moments = NaN(length(AR), 3);
% dev_moments = NaN(length(AR), 3);
% gam_moments = NaN(length(AR), 3);
% tw_moments = NaN(length(AR), 3);
% alpha = NaN(length(AR), 1);


%% Plot the eigenvalue distribution, CDF, spectrum, etc.
% hold on;
% for j = 1:length(AR)
% collection = [1, 3, 6, 9, 11, 13, 15, 17];
% omega = [0.0442, 0.0538, 0.0348, 0.0414, 0.0180, 0.0320, 0.0154, 0.0242];

%% what garch_eig_fourier_coef_T8e4.mat now contains
% [1:11, 13,15,17];

% collection = 1:10;
% omega = [0.0442, 0.0446, 0.0538, 0.0450, 0.0444, 0.0348, 0.0438, ...
%          0.0525, 0.0414, 0.0704];
collection = 1:length(AR);
n = 500;
% if exist('garch_eig_fourier_coef_T8e4.mat', 'file') == 2
%     load('garch_eig_fourier_coef_T8e4.mat');
% else
%     C = NaN(2*n+1, length(AR));
% end
% p = NaN(3, length(AR));
c = 1;
% maxeig = NaN(length(collection), fold);
% mineig = NaN(length(collection), fold);
% mag_params = NaN(2, length(collection));
% ph_params = NaN(3, length(collection));
Cii_params = NaN(4, 4);
% ev_params = NaN(9, 4);
% merge_point = NaN(9, 4);
N = 250;
%for q = [12, 16, 20, 24]
if 1
    
    %% Plot the maximum eigenvalue CDF
    % subplot(2, 2, c);
    load(sprintf('../data/GarchWishartT8e4Diag-%.3f.mat', 0), ...
         'Cii');
    Cii = reshape(Cii, prod(size(Cii)), 1);

    load(sprintf('../data/GarchWishartT8e4Eig-%.3f.mat', 0), ...
         'ev');
    ev = reshape(ev, prod(size(ev)), 1);
    Cii_params = stblfit(Cii, 'ecf');

    [y1, x1] = ecdf(Cii);
    y2 = stblcdf(x1, Cii_params(1), Cii_params(2), Cii_params(3), Cii_params(4));
    [y3, x3] = ecdf(ev);
    plot(log10(x3), log10(y3), log10(x1), log10(y1), log10(x1), ...
         log10(y2));
    legend('Eigenvalue', 'Cii', sprintf('S(%.2f, %.2f, %.2e, %.2e)', ...
                                        Cii_params(1), Cii_params(2), ...
                                        Cii_params(3), Cii_params(4))); 
    xlim([-2.3, -1.9]);
    ylim([-3.5, 0]);
    grid on
    % title(sprintf('N=250, T=%d', N*q));
    c = c + 1;
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

    %% Calculate the Fourier Transform
    % load(sprintf('GarchWishartT8e4Eig-%.3f.mat', AR(j)), ...
    %      'ev');
    % ev = reshape(ev, 1, prod(size(ev)));
    % cf = ecf(ev, [-n:n]);
    
    %% The phase
    % phase = angle(cf);
    % s = sign(phase);
    % s(s == 0) = 1;
    % s = (s(2:end) - s(1:end-1))/2;
    % I = find(s < 0) + 1;
    % I = [1, I, length(phase)+1];
    % L = length(I);
    % for m = L/2 + 1:L - 1
    %     phase(I(m):I(m+1)-1) = phase(I(m):I(m+1)-1) + 2*pi * (m-L/2);
    % end
    % for m = L/2:-1:2
    %     phase(I(m)-1:-1:I(m-1)) = phase(I(m)-1:-1:I(m-1)) - 2*pi * (L/2-m+1);
    % end
    % f = @(p, k) p(1) .* k + sign(k).*abs(p(2) .* k).^(p(3)).*tan(pi.*p(3)/2);
    % p = [2.743e-3 * tan(pi*AR(j)/2) + 8.557e-3, ...
    %      exp(0.3180*tan(AR(j)) - 7.165),...
    %      1.48];
    % P = lsqcurvefit(f, p, [-n:n], phase);
    % ph_params(:, c) = abs(P);

    %% The magnitude
    % mag = log(-log(abs(cf(n+2:end))));
    % x = log(1:n);
    % P = polyfit(x, mag, 1);
    % mag_params(:, c) = [P(1), exp(P(2)/P(1))];
    % %% plot
    % if c == 16
    %     figure;
    % end
    % if c <= 15
    %     subplot(3, 5, c);
    % else
    %     subplot(3, 4, c-15)
    % end
    % % plot([-n:n], phase, [-n:n], f(ph_params(:, c), [-n:n]));
    % plot(x, mag);
    % title(sprintf('\\phi = %.3f; P=%.2e %.2e', AR(j), abs(P)), ...
    %       'Interpreter', 'tex');
    % grid on
    % c = c + 1;

    %% fit Re c_n / Im c_n to dumped cosine / sine functions.
    % freal = @(p, x) exp(-p(1) - p(2) * abs(x)) .* cos(p(3) .* x);
    % q = lsqcurvefit(freal, [0.5, 4.5e-3, 6*pi/500], ...
    %                 [-n:n]', real(C(:, j)));
    % % fimag = @(p, x) exp(-p(1) - p(2) * abs(x)) .* sin(p(3) .* x);
    % % q = lsqcurvefit(fimag, [0.5, 4.5e-3, omega(j)], ...
    % %                 [-n:n]', imag(C(:, j)));
    % p(:, j) = q';
    % y = freal(q, [-n:n]');
    % subplot(4, 4, c);
    % plot(-n:n, real(C(:, j)), 'b', -n:n, y, 'g');
    % xlabel('n', 'Interpreter', 'tex');
    % ylabel('Im c_n', 'Interpreter', 'tex');
    % title(sprintf('\\phi=%.2f', AR(j)), 'Interpreter', 'tex');
    % grid on
    % c = c + 1;
    
    %% Compare the Fourier expansion with the empirical CDF.    
    % [y, x] = ecdf(ev);
    % x = reshape(x, length(x), 1);
    % K = exp(-i * (2*pi) * x * [-n:n] ./ (b-a)) - ...
    %     exp(-i * (2*pi) * (ones(length(x), 1).*a) * [-n:n] ./ (b-a));
    % K(:, n+1) = x - a;
    % F = (C(:, j) ./ [-n:n]') .* (i*(b-a)/(2*pi));
    % F(n+1) = C(n+1, j);
    % F = real(K * F);
    % F(F < 0) = 0;
    % plot(log(x), log(y), 'b', log(x), log(F), 'r');

    % xlabel('ln(\lambda)', 'Interpreter', 'tex');
    % ylabel('ln(F(\lambda))', 'Interpreter', 'tex');
    % title(sprintf('\\phi = %.2f', AR(collection(j))));

    %% Plot the phase of the Fourier coefficients
    % plot(-n:n, real(C(:, j)));
    % xlabel('n', 'Interpreter', 'tex');
    % ylabel('Re c_n', 'Interpreter', 'tex');
    % title(sprintf('\\phi=%.2f', AR(j)), 'Interpreter', 'tex');
    % grid on
    
    %% Plot the full eigen value spectrum
    % subplot(2, 2, c);
    % hold on

    % load(sprintf('GarchWishartN%dQ%dEig-%.3f.mat', N, q, 0), 'ev');
    % ev = reshape(ev, prod(size(ev)), 1);
    % [y1,x1] = ecdf(ev);
    
    % load(sprintf('GarchWishartN%dQ%dDiag-%.3f.mat', N, q, 0), 'Cii');
    % Cii = reshape(Cii, prod(size(Cii)), 1);
    % [y2,x2] = ecdf(Cii);
    
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

