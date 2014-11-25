clear all
close all

T = 8e+4;
N = 50;
fold = 2e3;

q = N/T;
AR = [0.825, 0.85, 0.875];
tau = -log(2)./log(AR);
% dist = struct('name', 'lorentzian', 'prmt', 1);
dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.15, 0.84], 'distr', ...
             'Gaussian', 'TailExponent', 2.9664);
% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.1, 0.725], 'distr', ...
%               struct('Name', 't', 'DoF', 3), 'TailExponent', 1.8187);
for j = AR
    clear ev Cij Cii;
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
    if (exist(sprintf('GarchWishartT8e4Eig-%.3f.mat', j), 'file') == 2)
        load(sprintf('GarchWishartT8e4Eig-%.3f.mat', j), 'ev');
        load(sprintf('GarchWishartT8e4Offdiag-%.3f.mat', j), 'Cij');
        load(sprintf('GarchWishartT8e4Diag-%.3f.mat', j), 'Cii');
        ev = [ev, ev1];
        Cij = [Cij, Cij1];
        Cii = [Cii, Cii1];
    else
        ev = ev1;
        Cij = Cij1;
        Cii = Cii1;
    end
    
    save(sprintf('GarchWishartT8e4Eig-%.3f.mat', j), 'ev');
    save(sprintf('GarchWishartT8e4Offdiag-%.3f.mat', j), 'Cij');
    save(sprintf('GarchWishartT8e4Diag-%.3f.mat', j), 'Cii');
    fprintf('Saved files %.3f.\n', j);
end
quit
