%% Generate a matrix of returns series with specfied AR. coefficients
function R = gen_ret_mtx (n, T, distr, coef, varargin)

% if strcmp(distr.name, 'garch')
%     model = garch('Constant', 2.3e-6, 'ARCH', 0.15, 'GARCH', ...
%                   0.84, 'Distribution', distr.prmt);
%     [V, R] = simulate(model, T, 'numPaths', n);
%     R = R';
%     return;
% end

% First find the auto-correlation coefficients
m = length(coef);
% R = NaN(n, T+m);

% x = NaN(n, T+m);
% x(:, 1:m) = zeros(n, m);
switch (distr.name)
  case {'Normal', 'Gaussian'}
    x = randn(n, T);
  case {'Lorentzian', 'Cauchy'}
    % A uniform distr. in (-1/2, 1/2) %
    x = rand(n, T) - 1/2;
    x = distr.prmt(1) .* tan(pi .* x);
    % y = randn(n, 2);
    % x = y(:, 1) ./ abs(y(:, 2));
  case 'Levy'
    alpha = distr.prmt(1);
    gam = distr.prmt(2);
    u = rand(n, T);
    v = rand(n, T);
    phi = pi * (v - 0.5);
    x = gam * (...
        -log(u) .* cos(phi) ./...
        cos((1 - alpha) * phi) ...
        ).^(1 - 1/alpha) .*...
        sin(alpha * phi) ./ cos(phi);
  case 'student t'
    dof = distr.prmt(1);
    y = randn(n, dof + 1);
    x = y(:, 1) .* sqrt(dof) ./ sqrt(sum(y(:, 2:end).^2, 2));
  case 'Garch1_1'
    a0 = distr.prmt(1);
    a1 = distr.prmt(2);
    b1 = distr.prmt(3);
    if a0 > 0 && (a1 + b1 < 1)
        h0 = a0/(1 - a1 - b1);
    else
        h0 = 1;
    end
    h = [ones(n, 1)*h0, NaN(n, T-1)];
    x = [randn(n, 1).*sqrt(h0), NaN(n, T-1)];
    for t = 2:T
        h(:, t) = a0 + a1.*x(:, t-1).^2 + b1.*h(:, t-1);
        x(:, t) = randn(n, 1).*sqrt(h(:, t));
    end
    
  case 'I-Garch1'
    a = distr.prmt;
    h0 = 1;
    h = [ones(n, 1)*h0, NaN(n, T-1)];
    x = [randn(n, 1).*sqrt(h0), NaN(n, T-1)];
    for t = 2:T
        h(:, t) = (1 - a).*x(:, t-1).^2 + a.*h(:, t-1);
        x(:, t) = randn(n, 1).*sqrt(h(:, t));
    end
end

if ~isempty(coef) && sum(coef) > 0
    R = [zeros(n, m), NaN(n, T)];
    for k = [1:T]
        R(:, k+m) = x(:, k) + R(:, [k-1:-1:k-m]+m) * diag(coef);
    end
    R = R(:, m+1:m+T);
else
    R = x;
end
