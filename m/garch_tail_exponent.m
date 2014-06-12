function k = garch_tail_exponent(dist, param, T);

if param(1) <= 0 || (param(2) + param(3) >= 1)
    k = NaN;
    return;
end
if strcmp(dist.name, 'Gaussian') == 1
    R = randn(1, T);
elseif strcmp(dist.name, 'Student t') == 1
    R = trnd(dist.DoF, 1, T);
end

A = param(2) .* R.^2 + param(3);
s = mean(log(A));

if s >= 0
    k = NaN;
    return;
end

a = 3;
while mean(A.^a) < 1
    a = a + 1;
end

f = @(k) mean(A.^(k/2)) - 1;
options = optimoptions('fsolve','Display','off');
k = fsolve(f, 2*a, options);
