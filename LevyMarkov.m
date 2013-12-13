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

for i = 1:length(phi)
    load(sprintf('MarkovEig-%d.mat', i-1), 'ev');
    maxeig = max(ev);
    [y, x] = ecdf(maxeig);
    plot(log(x), log(y), spec(i));
    hold on
end
hold off
