clear all
T = 1e3;
N = 50;
fold = 2e3;

q = N/T;

tau = 11:5:61;
AR = 2.^(-1./tau);
% AR = [[0:19]*0.05, 0.96:0.006:0.996, 0.9965, 0.997, 0.9975, 0.9980, ...
%       0.999];
% AR = [0:19]*0.05;
% AR = [0:2:18]*0.05;
for j = 1:length(AR)
    if (exist(sprintf('GaussianWishartOffDiagLinearTau-%d.mat', tau(j)), 'file') == 2)
        continue;
    end
    % ev = NaN(N, fold);
    offdiag = NaN(N*(N-1), fold);
    for i = 1:fold
        A = randn(N, T);
        A1 = zeros(N, T+1);
        if AR(j) ~= 0
            for k = 2:T+1
                A1(:, k) = A1(:, k-1)*AR(j) + A(:, k-1);
            end
            A1 = A1(:, 2:end);
        else
            A1 = A;
        end
        % S = std(A1');
        % A1 = diag(1./S)*A1;
        C = A1*A1'./T;
        % ev(:, i) = eig(C);
        offdiag(:, i) = C(~eye(N));
    end
    save(sprintf('GaussianWishartOffDiagLinearTau-%d.mat', tau(j)), 'offdiag');
    fprintf('Saved file %d.\n', tau(j));
end

mu = NaN(length(tau), 1);
v = NaN(length(tau), 1);
for j = 1:length(tau)
    load(sprintf('GaussianWishartOffDiagLinearTau-%d.mat', tau(j)), ...
         'offdiag');
    offdiag = reshape(offdiag, prod(size(offdiag)), 1);
    mu(j) = mean(offdiag);
    v(j) = std(offdiag);
end
