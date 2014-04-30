clear all
close all

N = 500;
q = 12;
T = N * q;
fold = 2000;


AR = 0:0.1:0.9;
% AR = 0;
tau = -log(2)./log(AR);
dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.90, 0.09], 'distr', ...
             struct('Name', 'Gaussian'), 'TailExponent', 2.0329);
% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.15, 0.84], 'distr', ...
%              'Gaussian', 'TailExponent', 2.9664);
% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.1, 0.8], 'distr', ...
%               struct('Name', 't', 'DoF', 3), 'TailExponent', 0.7744);
for j = AR
    ev1 = NaN(N, fold);
    Cij1 = NaN(N*(N-1)/2, fold);
    Cii1 = NaN(N, fold);
    for i = 1:fold
        R = gen_ret_mtx(N, T, dist, j);
        if strcmp(dist.name, 'Cauchy') == 1
            R = R ./ dist.prmt ./ T;
        elseif strcmp(dist.name, 'Garch1_1') == 1
            R = R ./ T^(1/dist.TailExponent);
        end
        C = R*R';
        Cij1(:, i) = C(logical(triu(ones(N), 1)));
        Cii1(:, i) = C(logical(eye(N)));
        ev1(:, i) = eig(C);
    end
    if (exist(sprintf('GarchWishart_0.90_0.09N%dQ%dEig-%.3f.mat', ...
                      N, q, j), 'file') == 2)
        load(sprintf('GarchWishart_0.90_0.09N%dQ%dEig-%.3f.mat', ...
                     N, q, j), 'ev');
        load(sprintf(['GarchWishart_0.90_0.09N%dQ%dOffdiag-' ...
        '%.3f.mat'], N, q, j), 'Cij');
        load(sprintf('GarchWishart_0.90_0.09N%dQ%dDiag-%.3f.mat', ...
                     N, q, j), 'Cii');
        ev = [ev, ev1];
        Cij = [Cij, Cij1];
        Cii = [Cii, Cii1];
    else
        ev = ev1;
        Cij = Cij1;
        Cii = Cii1;
    end
    
    save(sprintf('GarchWishart_0.90_0.09N%dQ%dEig-%.3f.mat', N, q, j), 'ev');
    save(sprintf('GarchWishart_0.90_0.09N%dQ%dOffdiag-%.3f.mat', N, q, j), 'Cij');
    save(sprintf('GarchWishart_0.90_0.09N%dQ%dDiag-%.3f.mat', N, q, j), 'Cii');
    fprintf('Saved files %.3f.\n', j);
end
quit
