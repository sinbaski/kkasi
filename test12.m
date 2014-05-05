clear all
close all

fold = 4e+3;
% fold = 1;

%AR = 0:0.1:0.9;
AR = [0.2:0.2:0.8];
% AR = 0;
tau = -log(2)./log(AR);
% dist = struct('name', 'Gaussian', 'TailExponent', 2);
% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.15, 0.84], 'distr', ...
%               struct('Name', 'Gaussian'), 'TailExponent', 2.9664);
dist = struct('name', 'Garch1_1', 'prmt', [0.01, 0.15, 0.84], 'distr', ...
              struct('Name', 'Gaussian'), 'TailExponent', 2.9664);
% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.1, 0.725], 'distr', ...
%               struct('Name', 't', 'DoF', 3), 'TailExponent',
%               1.8187);
%for q = [12, 16, 20, 24]
N = 2;
% theta = pi/4;
rho = 0;
C = eye(N);
C(~logical(eye(N))) = rho;
[U, D] = eig(C);
W = U*sqrt(D);

for q = 1500
    T = N * q;
    for j = AR
        % if (exist(sprintf('GarchCovN%dQ%dDiag-%.3f.mat', N, ...
        %                   T, j), 'file') == 2)
        %     continue;
        % end
        ev1 = NaN(N, fold);
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
            % R = diag(1./std(R')) * R;
            if rho == 0
                % C1(:, :, i) = R*R' ./ T^(2/dist.TailExponent);
                C1(:, :, i) = R*R' ./ T;
            else
                Z = W * R;
                % C1(:, :, i) = Z*Z' ./ T^(2/dist.TailExponent);
                C1(:, :, i) = Z*Z' ./ T;
            end
            [V, D] = eig(C1(:, :, i));
            ev1(:, i) = D(logical(eye(N)));
            X1(:, :, i) = V;
        end
        if (exist(sprintf('GarchCovN%dQ%dRho%.3fEig-%.3f.mat', ...
                     N, q, rho, j), 'file') == 2)
            % load(sprintf('GarchCovN%dQ%dRho%.3fEig-%.3f.mat', ...
            %              N, q, rho, j), 'ev');
            % load(sprintf('GarchCovN%dQ%dRho%.3fX-%.3f.mat', ...
            %              N, q, rho, j), 'X');
            % load(sprintf('GarchCovN%dQ%dRho%.3fC-%.3f.mat', ...
            %              N, q, rho, j), 'C');
            % ev = [ev, ev1];
            % X = cat(3, X, X1);
            % C = cat(3, C, C1);
        else
            % ev = ev1;
            % X = X1;
            % C = C1;
        end
        ev = ev1;
        X = X1;
        C = C1;

        save(sprintf('GarchCovN%dQ%dRho%.3fEig-%.3f.mat', ...
                     N, q, rho, j), 'ev');
        save(sprintf('GarchCovN%dQ%dRho%.3fX-%.3f.mat', ...
                     N, q, rho, j), 'X');
        save(sprintf('GarchCovN%dQ%dRho%.3fC-%.3f.mat', ...
                     N, q, rho, j), 'C');
        fprintf('Saved files %.3f.\n', j);
    end
end
% quit
