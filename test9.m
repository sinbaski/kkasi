clear all
% f = @(x) A*x.^3 + B*x + 1;

% x = -1:1e-2:1;
% y = erfc(x);

% plot(x, y);
T = 1e3;
N = 50;
fold = 2e3;

% AR = [[0:19]*0.05, 0.96:0.006:0.996, 0.9965, 0.997, 0.9975, 0.9980, ...
%       0.999];
AR = [0:19]*0.05;

q = N/T;
d = 20;
% phi = 0.05 * (d - 1);
phi = AR(d);
% load(sprintf('GaussianWishartEig005-%d.mat', d), 'ev');
load(sprintf('GaussianWishartOffDiagLinearTau-%d.mat', d), 'ev');
for i = 1:fold
    A = randn(N, T);
    A1 = zeros(N, T+1);
    if phi ~= 0
        for k = 2:T+1
            A1(:, k) = A1(:, k-1)*phi + A(:, k-1);
        end
        A1 = A1(:, 2:end);
    else
        A1 = A;
    end
    % S = std(A1');
    % A1 = diag(1./S)*A1;
    C = A1*A1'./T;
    ev1(:, i) = eig(C);
end
ev = [ev, ev1];
save(sprintf('GaussianWishartOffDiagLinearTau-%d.mat', d), 'ev');
