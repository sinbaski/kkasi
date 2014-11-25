clear all
close all

N = 50;
T = 1e3;
fold = 2000;
alpha = 1;
gam = 1;
% dist = struct('name', 'Levy', 'prmt', [alpha, gam]);
dist = struct('name', 'lorentzian', 'prmt', 1);
phi = [0:5] * 0.15;
spec = ['b', 'c', 'g', 'm', 'r', 'k'];
for i = 1:length(phi)
    if (exist(sprintf('./MarkovEig-%d.mat', i-1), 'file') == 2)
        continue;
    end
    ev = NaN(N, fold);
    Cii = NaN(N, fold);
    for j = 1:fold
        R = gen_ret_mtx(N, T, dist, [phi(i)]);
        R = R./T;
        % s = diag(1./std(R'));
        % M = s*R;
        C = R*R';
        ev(:, j) = eig(C);
        for k = 1:N
            Cii(k, j) = C(k, k);
        end
    end
    save(sprintf('MarkovEig-%d.mat', i-1), 'ev', 'Cii');
end

figure;
for i = 1:length(phi)
    load(sprintf('MarkovEig-%d.mat', i-1), 'ev');
    % maxeig = max(ev);
    %% CDF
    % [y, x] = ecdf(maxeig);
    % plot(x, log(y), spec(i));

    %% PDF
    % lag = 16;
    % x = linspace(0.12, 0.24, 160);
    % y = hist(maxeig, x) ./ length(maxeig) ./ (x(2) - x(1));
    % x1 = tsmovavg(x, 's', lag);
    % x1 = x1(lag:end);
    % y1 = tsmovavg(y, 's', lag);
    % y1 = y1(lag:end);
    % plot(x1, y1, spec(i));
    
    %% spectrum
    lag = 16;
    ev = reshape(ev, N*fold, 1);
    % x = linspace(min(ev), max(ev), 200);
    x = linspace(min(ev), 5, 450);
    % x = [x(1) - (x(2)-x(1))*[lag-1:-1:1], x];
    y = hist(ev, x) ./ (N*fold) ./ (x(2)-x(1));
    y = y(1:end-1);
    x = x(1:end-1);
    x = tsmovavg(x, 's', lag);
    y = tsmovavg(y, 's', lag);
    x = x(lag:end);
    y = y(lag:end);
    plot(x, y, spec(i));
    hold on
end
hold off
grid on
title(sprintf('PDF of eigen values. a_t ~ Cauchy(1). N/T=%d/%d', ...
              N, T));
ylabel('$f(\lambda)$', 'Interpreter', 'Latex', 'Fontsize', 14);
xlabel('$\lambda$',  'Interpreter', 'Latex', 'Fontsize', 14);

legends = cell(1, length(phi));
for i = 1:length(legends)
    legends{i} = sprintf('r_t = %.2f r_{t-1} + a_t', phi(i));
end
legend(legends{1}, legends{2}, legends{3}, legends{4}, legends{5}, ...
       legends{6}, 'Location', 'Northeast');

%% Plot the diagonal elements distribution
figure;
for i = 1:length(phi)
    load(sprintf('MarkovEig-%d.mat', i-1), 'Cii');
    lag = 16;
    Cii = reshape(Cii, N*fold, 1);
    x = linspace(min(Cii), 5, 450);
    % x = [x(1) - (x(2)-x(1))*[lag-1:-1:1], x];
    y = hist(Cii, x) ./ (N*fold) ./ (x(2)-x(1));
    y = y(1:end-1);
    x = x(1:end-1);
    x = tsmovavg(x, 's', lag);
    y = tsmovavg(y, 's', lag);
    x = x(lag:end);
    y = y(lag:end);
    plot(x, y, spec(i));
    hold on
end
hold off
grid on
title(sprintf('PDF of diagonal elements. a_t ~ Cauchy(1). N/T=%d/%d', ...
              N, T));
ylabel('$P(\phi; C_{ii})$', 'Interpreter', 'Latex', 'Fontsize', 14);
xlabel('$C_{ii}$',  'Interpreter', 'Latex', 'Fontsize', 14);

legends = cell(1, length(phi));
for i = 1:length(legends)
    legends{i} = sprintf('r_t = %.2f r_{t-1} + a_t', phi(i));
end
legend(legends{1}, legends{2}, legends{3}, legends{4}, legends{5}, ...
       legends{6}, 'Location', 'Northeast');
