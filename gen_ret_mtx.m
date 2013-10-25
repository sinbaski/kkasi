% Generate a matrix of returns series.
function R = gen_ret_mtx (n, T, distr, auto_corr)

if strcmp(distr.name, 'garch')
    model = garch('Constant', 2.3e-6, 'ARCH', 0.15, 'GARCH', ...
                  0.84, 'Distribution', distr.prmt);
    [V, R] = simulate(model, T, 'numPaths', n);
    R = R';
    return;
end

% First find the auto-correlation coefficients
coef = [];
m = 1;
x = auto_corr(1);
while abs(x) > 0.05
    coef(m) = x;
    m = m + 1;
    x = auto_corr(m);
end

coef = flipud(coef');

R = zeros(n, T);
for col = 1:T
    switch (distr.name)
      case 'normal'
        x = randn(n, 1);
      case 'lorentzian'
        % A uniform distr. in (-1/2, 1/2) %
        x = rand(n, 1) - 1/2;
        x = distr.prmt(1) .* tan(pi .* x);
        % y = randn(n, 2);
        % x = y(:, 1) ./ abs(y(:, 2));
      case 'Levy'
        alpha = distr.prmt(1);
        gamma = distr.prmt(2);
        u = rand(n, 1);
        v = rand(n, 1);
        phi = pi * (v - 0.5);
        x = gamma * (...
            -log(u) .* cos(phi) ./...
            cos((1 - alpha) * phi) ...
            ).^(1 - 1/alpha) .*...
        sin(alpha * phi) ./ cos(phi);
      case 'student t'
        dof = distr.prmt(1);
        y = randn(n, dof + 1);
        x = y(:, 1) .* sqrt(dof) ./ sqrt(sum(y(:, 2:end).^2, 2));
    end
    if col == 1 || isempty(coef)
        R(:, col) = x;
    else
        m = min(col - 1, length(coef));
        w = coef(end - m + 1 : end);
        R(:, col) = R(:, col - m : col - 1) * w + ...
            (1 - sum(w)) .* x;
    end
end

