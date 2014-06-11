clear all
close all

fold = 4e+3;
% fold = 1;

%AR = 0:0.1:0.9;
AR = [0:0.2:0.8];
% AR = 0;
tau = -log(2)./log(AR);
dist = struct('name', 'Gaussian', 'TailExponent', 2);
% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.15, 0.84], 'distr', ...
%               struct('Name', 'Gaussian'), 'TailExponent', 2.9664);
% dist = struct('name', 'Garch1_1', 'prmt', [0.01, 0.15, 0.84], 'distr', ...
%               struct('Name', 'Gaussian'), 'TailExponent', 2.9664);
% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.1, 0.725], 'distr', ...
%               struct('Name', 't', 'DoF', 3), 'TailExponent',
%               1.8187);
%for q = [12, 16, 20, 24]
N = 4;
% theta = pi/4;
rho = 0;
C = eye(N);
C(~logical(eye(N))) = rho;
[U, D] = eig(C);
W = U*sqrt(D);

for q = 1500
    T = N * q;
    for j = AR
        % if (exist(sprintf('GarchWishartN%dQ%dDiag-%.3f.mat', N, ...
        %                   T, j), 'file') == 2)
        %     continue;
        % end
        ev1 = NaN(N, fold);
        % sig1 = NaN(N, fold);
        X1 = NaN(N, N, fold);
        C1 = NaN(N, N, fold);
        for i = 1:fold
            R = gen_ret_mtx(N, T, dist, j);
            % if strcmp(dist.name, 'Cauchy') == 1
            %     R = R ./ dist.prmt ./ T;
            % elseif strcmp(dist.name, 'Garch1_1') == 1
            %     % R = R ./ T^(1/dist.TailExponent);
            %     ;
            % end
            if rho == 0
                C1(:, :, i) = R*R' ./ T;
            else
                Z = W * R;
                C1(:, :, i) = Z*Z' ./ T;
            end
            [V, D] = eig(C1(:, :, i));
            ev1(:, i) = D(logical(eye(N)));
            X1(:, :, i) = V;
        end
        if (exist(sprintf('GaussianWishartN%dQ%dRho%.3fEig-%.3f.mat', ...
                     N, q, rho, j), 'file') == 2)
            load(sprintf('GaussianWishartN%dQ%dRho%.3fEig-%.3f.mat', ...
                         N, q, rho, j), 'ev');
            % load(sprintf('GaussianWishartN%dQ%dRho%.3fSig-%.3f.mat', ...
            %              N, q, rho, j), 'sig');
            load(sprintf('GaussianWishartN%dQ%dRho%.3fX-%.3f.mat', ...
                         N, q, rho, j), 'X');
            load(sprintf('GaussianWishartN%dQ%dRho%.3fC-%.3f.mat', ...
                         N, q, rho, j), 'C');
            ev = [ev, ev1];
            % sig = [sig, sig1];
            X = cat(3, X, X1);
            C = cat(3, C, C1);
        else
            ev = ev1;
            % sig = sig1;
            X = X1;
            C = C1;
        end

        save(sprintf('GaussianWishartN%dQ%dRho%.3fEig-%.3f.mat', ...
                     N, q, rho, j), 'ev');
        % save(sprintf('GaussianWishartN%dQ%dRho%.3fSig-%.3f.mat', ...
        %              N, q, rho, j), 'sig');
        save(sprintf('GaussianWishartN%dQ%dRho%.3fX-%.3f.mat', ...
                     N, q, rho, j), 'X');
        save(sprintf('GaussianWishartN%dQ%dRho%.3fC-%.3f.mat', ...
                     N, q, rho, j), 'C');
        fprintf('Saved files %.3f.\n', j);
    end
end
quit
