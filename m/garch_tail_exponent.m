function k = garch_tail_exponent(dist, param, T);
% function k = garch_tail_exponent(param);

a1 = param(1);
b1 = param(2);

if strcmp(dist.name, 'Gaussian') == 1
    R = randn(1, T);
elseif strcmp(dist.name, 'Student t') == 1
    R = trnd(dist.DoF, 1, T);
end
if param(1) <= 0
    k = NaN;
    return;
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
k = fsolve(f, 2*a);
