clear all
close all

T = 1e3;
N = 50;
fold = 2e3;

q = N/T;
ev = NaN(1, N*fold);
% mp_max = (1 + sqrt(q))^2;

%% eigen values generated from the volvo stochastic log-volatility
%% model.

spec = cellstr(['b  '; 'g  '; 'r  ';]);
% spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
%                 'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
%                 'b-.'; 'c-.'; 'g-.'; 'm-.'; 'r-.'; 'k-.';...
%                 ]);
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
%      0.999];
AR = [0:19]*0.05;
tau = -log(2)./log(AR);
% AR = [0:2:18]*0.05;

% for j = 1:length(AR)
%     if (exist(sprintf('GaussianWishartEig-%.3f.mat', AR(j)), 'file') == 2)
%         continue;
%     end
%     ev = NaN(N, fold);
%     for i = 1:fold
%         A = randn(N, T);
%         A1 = zeros(N, T+1);
%         if AR(j) ~= 0
%             for k = 2:T+1
%                 A1(:, k) = A1(:, k-1)*AR(j) + A(:, k-1);
%             end
%             A1 = A1(:, 2:end);
%         else
%             A1 = A;
%         end
%         C = A1*A1'./T;
%         ev(:, i) = eig(C);
%     end
%     save(sprintf('GaussianWishartEig005-%d.mat', j), 'ev');
%     fprintf('Saved file %d.\n', j);
% end


% [y, x] = ecdf(maxeig);
% p = RealWishartCDF(x(1:10:end), N, T);

% plot(x, y, x(1:10:end), p);

%% Approximation of the Tracy-Widom with a Gamma dist.
a1 = -1/2;
a2 = a1;
n = N;
p = T;
mu = (sqrt(n+a1) + sqrt(p+a2))^2;
sigma = sqrt(mu) * (1/sqrt(n+a1) + 1/sqrt(p+a2))^(1/3);

% hold on;
% avgs = NaN(1, length(AR));

sample_moments = NaN(length(AR), 3);
% dev_moments = NaN(length(AR), 3);
gam_moments = NaN(length(AR), 3);
tw_moments = NaN(length(AR), 3);
alpha = NaN(length(AR), 1);
k = NaN(length(AR), 1);
theta = NaN(length(AR), 1);


%% Plot the eigenvalue distribution, CDF, spectrum, etc.
hold on;
mineig = NaN(1, length(AR));
c = 1;
for j = [0, 0.4, 0.8]
% for j = 1:length(AR)
    load(sprintf('../matfys/data/GaussianWishartEig-%.3f.mat', j), 'ev');
    % mineig(j) = mean(min(ev));

    % maxeig = max(ev);
    % sample_moments(j, :) = [mean(maxeig), std(maxeig), ...
    %                     skewness(maxeig)];
    
    % maxeig1 = (maxeig - mu)/sigma;
    % tw_moments(j, :) = [mean(maxeig1), var(maxeig1), ...
    %                     skewness(maxeig1)];
    % moments = tw_moments(j, :);

    % k(j) = 4/(moments(3)^2);
    % theta(j) = sqrt(moments(2)) * moments(3)/2;
    % alpha(j) = k(j)*theta(j) - moments(1);
    % maxeig1 = maxeig1 + alpha(j);
    % gam_moments(j, :) = [mean(maxeig1), var(maxeig1), ...
    %                     skewness(maxeig1)];
    
    %% Plot the full eigen value spectrum
    ev = reshape(ev, prod(size(ev)), 1);
    lag = 4;
    x = linspace(min(ev)*0.9, max(ev)*1.05, 240);
    y = hist(ev, x) / length(ev) / (x(2) - x(1));
    x = tsmovavg(x, 's', lag);
    x = x(lag:end);
    y = tsmovavg(y, 's', lag);
    y = y(lag:end);
    if c == 0
        y1 = MarcenkoPasturPDF(x, [N/T, 1]);
        stairs(x, y1);
    end
    plot(x, y, spec{c}, 'LineWidth', 2);
    c = c + 1;
    
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
    % y1 = cdf('Gamma', x, k(j), theta(j));
    % subplot(2, 3, c);
    % plot(log(x), log(y), 'b', log(x), log(y1), 'r');
    % title(sprintf('\\phi=%.3f. \\Gamma(%.1f, %.2e)', AR(j), ...
    %               k(j), theta(j)));
    % % legend(sprintf('\\phi=%.3f', AR(j)),...
    % %        sprintf('\\Gamma(%.1f, %.2e)', k(j), theta(j)), 'Location', ...
    % %        'Southeast');
    % grid on
    % xlabel('ln((\lambda_1 - \mu_{NT})/\sigma_{NT})');
    % ylabel('ln(F(\lambda_1))');
    % c = c + 1;
    %% plot CDF of maxeig - (1 + sqrt(q))^2
    % dev = maxeig - (1 + sqrt(q))^2 / (1 - AR(j)^2);
    % dev_moments(j, :) = [mean(dev), std(dev), skewness(dev)];
    % [y, x] = ecdf(dev);
    % y1 = cdf('Normal', x, mean(dev), std(dev));
    % plot(x, y, x, y1, spec{j});

end
hold off
grid on

% x = tan(pi.*AR(1:20)'./2);
% y1 = theta(1:20);
% y2 = k(1:20).*theta(1:20);
% y3 = k(1:20).*theta(1:20).^2;

% p1 = polyfit(x, y1, 2);
% r1 = NaN(20, 1);
% for n = 1:length(r1)
%     z = roots([p1(1:2), p1(3) - y1(n)]);
%     r1(n) = z(imag(z) == 0 & real(z) > 0);
% end

% p2 = polyfit(x, y2, 2);
% r2 = (-p2(2) + sqrt(p2(2)^2 - 4.*p2(1).*(p2(3) - y2)))...
%      ./(2 * p2(1));


% p3 = polyfit(x, y3, 3);
% r3 = NaN(20, 1);
% for n = 1:length(r3)
%     z = roots([p3(1:3), p3(4) - y3(n)]);
%     r3(n) = z(imag(z) == 0);
% end

% subplot(2, 2, 1);
% plot(AR, k.*theta, 'b+-');
% grid on
% xlabel('Auto-correlation \phi');
% ylabel('mean of \lambda_1');

% subplot(2, 2, 2);
% plot(AR, k.*theta.^2, 'b+-');
% grid on
% xlabel('Auto-correlation \phi');
% ylabel('variance of \lambda_1');

% subplot(2, 2, 3);
% plot(AR, 2./sqrt(k), 'b+-');
% grid on
% xlabel('Auto-correlation \phi');
% ylabel('variance of \lambda_1');

% subplot(2, 2, 4);
% plot(AR, alpha, 'b+-');
% grid on
% xlabel('Auto-correlation \phi');
% ylabel('Relocation Parameter \alpha');



% grid on
% hold off;
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

% legend('$\tau$=0; theoretical', '$\tau$=0', '$\tau$=0.30', '$\tau$=0.43', ...
%        '$\tau$=0.58', '$\tau$=0.75', '$\tau$=1.0', '$\tau$=1.36', '$\tau$=1.94', ...
%        '$\tau$=3.11', 'Location', 'Northeast');

% % legend('$\tau=11$', '$\tau=21$', '$\tau=31$', '$\tau=41$');
% h = legend;
% set(h, 'interpreter', 'latex', 'fontsize', 14);
% xlabel('$\lambda_1$', 'Interpreter', 'LaTex');
% ylabel('$f_{\tau}(\lambda_1)$', 'Interpreter', 'latex', 'Fontsize', 14);
