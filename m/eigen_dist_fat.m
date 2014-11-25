clear all
close all

T = 1e3;
N = 50;
fold = 2e3;

q = N/T;

%% eigen values generated from the volvo stochastic log-volatility
%% model.

spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
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

% tau = 11:5:61;
% AR = 2.^(-1./tau);
% tau = [0, tau];
% AR = [0, AR];
% AR = [[0:19]*0.05, 0.96:0.006:0.996, 0.9965, 0.997, 0.9975, 0.9980, ...
%       0.999];
AR = [[0:9]*0.1, 0.95, 0.99];
tau = -log(2)./log(AR);
% AR = [0:2:18]*0.05;

for j = AR
    if (exist(sprintf('FatWishartDiag-%.3f.mat', j), 'file') == 2)
        continue;
    end
    ev = NaN(N, fold);
    Cij = NaN(N*(N-1)/2, fold);
    Cii = NaN(N, fold);
    for i = 1:fold
        A = randn(N, T, 2);
        R = A(:, :, 1) .* exp(A(:, :, 2));
        if j ~= 0
            for t = 2:T
                R(:, t) = R(:, t-1)*j + R(:, t);
            end
        end
        R = diag(1./std(R')) * R ./ sqrt(T);
        C = R*R';
        Cij(:, i) = C(logical(triu(ones(N), 1)));
        Cii(:, i) = C(logical(eye(N)));
        ev(:, i) = eig(C);
    end
    save(sprintf('FatWishartEig-%.3f.mat', j), 'ev');
    save(sprintf('FatWishartOffdiag-%.3f.mat', j), 'Cij');
    save(sprintf('FatWishartDiag-%.3f.mat', j), 'Cii');
    fprintf('Saved file %.3f.\n', j);
end


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
c = 1;
n = 500;
collection = 1:length(AR);

if exist('fat_eig_fourier_coef.mat', 'file') == 2
    load('fat_eig_fourier_coef.mat');
else
    C = NaN(2*n+1, length(collection));
end
p = NaN(3, length(collection));

for j = collection
    %% Calculate the Fourier coefficients
    % load(sprintf('FatWishartEig-%.3f.mat', AR(j)), 'ev');
    % ev = reshape(ev, 1, prod(size(ev)));
    % a = min(ev);
    % b = max(ev);
    % C(:, j) = fourierExpandComplex(ev, n, b - a);

    subplot(3, 4, c);
    c = c + 1;

    %% fit Re c_n / Im c_n to dumped cosine / sine functions.
    % freal = @(p, x) exp(-p(1) - p(2) * abs(x).^(1/2)) .* cos(p(3) .* x);
    % q = lsqcurvefit(freal, [1, 1, pi], [-n:n]', real(C(:, j)));
    % % fimag = @(p, x) exp(-p(1) - p(2) * abs(x)) .* sin(p(3) .* x);
    % % q = lsqcurvefit(fimag, [0.5, 4.5e-3, omega(j)], ...
    % %                 [-n:n]', imag(C(:, j)));
    % p(:, j) = q';
    % y = freal(q, [-n:n]');
    % plot(-25:25, real(C([-25:25]+501, j)), 'b', [-25:25], y([-25:25]+501), 'g');
    % xlabel('n', 'Interpreter', 'tex');
    % ylabel('Re c_n', 'Interpreter', 'tex');
    % title(sprintf('\\phi=%.2f', AR(j)), 'Interpreter', 'tex');
    % grid on
 
    %% Plot the Fourier coefficients
    x = -n:n;
    fnorm = @(p, u) p(1) ./ (abs(u) + p(2)).^(p(3));
    q = lsqcurvefit(fnorm, [1, 1, 1], x', abs(C(x+501, j)));
    p(:, j) = q';
    plot(x, abs(C(x+501, j)), 'b', x, fnorm(q, x));
    xlabel('n', 'Interpreter', 'tex');
    ylabel('|c_n|', 'Interpreter', 'tex');
    title(sprintf('\\phi=%.2f', AR(collection(j))), 'Interpreter', 'tex');
    grid on

    %% Plot the phase of the Fourier coefficients
    % plot(-n:n, real(C(:, j)), 'b');
    % xlabel('n', 'Interpreter', 'tex');
    % ylabel('ph c_n', 'Interpreter', 'tex');
    % title(sprintf('\\phi=%.2f', AR(collection(j))), 'Interpreter', 'tex');
    % grid on
    % c = c + 1;
    
    
    %% Plot the full eigen value spectrum
    % ev = reshape(ev, prod(size(ev)), 1);
    % lag = 4;
    % [x, y] = epdf(ev, lag, 0, 6, 400, spec{c});
    % if AR(j) == 0
    %     y1 = MarcenkoPasturPDF([N/T, 1], x);
    %     stairs(x, y1, 'r-.', 'LineWidth', 1);
    % end
    % c = c + 1;
    
    %% Plot the PDF
    % lag = 16;
    % x = linspace(min(maxeig)*0.9, max(maxeig)*1.1, 120);
    % y = hist(maxeig, x) / length(maxeig) / (x(2) - x(1));
    % x = tsmovavg(x, 's', lag);
    % x = x(lag:end);
    % y = tsmovavg(y, 's', lag);
    % y = y(lag:end);
    % % y1 = pdf('Gamma', x, k, theta);
    % % plot(x, y, spec{(j-1)/4+1});
    % plot(x, y, spec{j});

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

end
if exist('fat_eig_fourier_coef.mat', 'file') ~= 2
    save('fat_eig_fourier_coef.mat', 'C');
end
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

% legend('$\phi=0$', 'MP law (q=0.05, $\sigma$=1)', '$\phi=0.5$', ...
%        '$\phi=0.8$', '$\phi=0.95$', 'Location', 'Northeast');
% h = legend;
% set(h, 'interpreter', 'latex', 'fontsize', 14);
% xlabel('$\lambda$', 'Interpreter', 'LaTex');
% ylabel('$f_{\phi}(\lambda)$', 'Interpreter', 'latex', 'Fontsize', 14);
