clear all
close all

N = 50;
T = 1e3;
fold = 2000;
alpha = 3/2;
gam = 1;
dist = struct('name', 'Levy', 'prmt', [alpha, gam]);
phi = [0:5] * 0.15;
spec = ['b', 'c', 'g', 'm', 'r', 'k'];
for i = 1:length(phi)
    if (exist(sprintf('./MarkovEig-%d.mat', i-1), 'file') == 2)
        continue;
    end
    ev = NaN(N, fold);
    for j = 1:fold
        R = gen_ret_mtx(N, T, dist, [phi(i)]);
        s = diag(1./std(R'));
        M = s * R;
        C = T.^(-2/alpha) * M*M';
        ev(:, j) = eig(C);
    end
    save(sprintf('MarkovEig-%d.mat', i-1), 'ev');
end

figure;
for i = 1:length(phi)
    load(sprintf('MarkovEig-%d.mat', i-1), 'ev');
    maxeig = max(ev);
    %% CDF
    % [y, x] = ecdf(maxeig);
    % plot(x, log(y), spec(i));

    %% PDF
    lag = 16;
    x = linspace(0.12, 0.24, 160);
    y = hist(maxeig, x) ./ length(maxeig) ./ (x(2) - x(1));
    x1 = tsmovavg(x, 's', lag);
    x1 = x1(lag:end);
    y1 = tsmovavg(y, 's', lag);
    y1 = y1(lag:end);
    plot(x1, y1, spec(i));
    hold on
end
hold off
grid on
title(sprintf('PDF of Maximum eigen value. a_t ~ Levy(%.2f, %.2f). N/T=%d/%d', ...
              alpha, gam, N, T));
ylabel('f(\lambda_{max})');
xlabel('\lambda_{max}');

legends = cell(1, length(phi));
for i = 1:length(legends)
    legends{i} = sprintf('r_t = %.2f r_{t-1} + a_t', phi(i));
end
legend(legends{1}, legends{2}, legends{3}, legends{4}, legends{5}, ...
       legends{6}, 'Location', ...
       'Southeast');
