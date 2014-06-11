function k = garch_tail_exponent(dist, a1, b1);

T = 4e+6;
if strcmp(dist.name, 'Gaussian') == 1
    R = randn(1, T);
elseif strcmp(dist.name, 'Student t') == 1
    R = trnd(dist.DoF, 1, T);
end
A = a1 .* R.^2 + b1;
f = @(k) mean(A.^(k/2)) - 1;
k = fsolve(f, 3);
