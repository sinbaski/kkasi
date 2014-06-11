clear all
close all

T = 1e3;
N = 50;
fold = 3e5;

q = N/T;
spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'y  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'y--'; 'k--';...
                ]);
AR = 0.5;
% ev = NaN(N, fold);
% for i = 1:fold
%     A = randn(N, T);
%     A1 = zeros(N, T+1);
%     if AR ~= 0
%         for k = 2:T+1
%             A1(:, k) = A1(:, k-1)*AR + A(:, k-1);
%         end
%         A1 = A1(:, 2:end);
%     else
%         A1 = A;
%     end
%     % S = std(A1');
%     % A1 = diag(1./S)*A1;
%     C = A1*A1'./T;
%     ev(:, i) = eig(C);
% end
% save(sprintf('GaussianWishartEig-6-extra5.mat'), 'ev');
ev1 = [];
for i = 1:5
    load(sprintf('GaussianWishartEig-6-extra%d.mat', i), 'ev');
    ev1 = [ev1, ev];
end
ev = ev1;

a1 = -1/2;
a2 = a1;
n = N;
p = T;
mu = (sqrt(n+a1) + sqrt(p+a2))^2;
sigma = sqrt(mu) * (1/sqrt(n+a1) + 1/sqrt(p+a2))^(1/3);
sample_moments = NaN(length(AR), 3);
dev_moments = NaN(length(AR), 3);
gam_moments = NaN(length(AR), 3);

maxeig = max(ev);
sample_moments(1, :) = [mean(maxeig), std(maxeig), ...
                    skewness(maxeig)];

maxeig1 = (maxeig - mu)/sigma;
moments = [mean(maxeig1), var(maxeig1), skewness(maxeig1)];

k = 4/(moments(3)^2);
theta = sqrt(moments(2)) * moments(3)/2;
alpha = k*theta - moments(1);
maxeig1 = maxeig1 + alpha;
gam_moments(1, :) = [mean(maxeig1), std(maxeig1), ...
                    skewness(maxeig1)];
idx = 501:800;
plot_pdf(maxeig1, 800, @(x) pdf('Gamma', x, k, theta), 15, idx, 1, '', '');
