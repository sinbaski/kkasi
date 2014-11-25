function F = LognormalEndsTailFun(lam1, lam_max, v, q, a, C)
if isempty(C)
    C = NaN(1,3);
    a = LognormalEndsConstA(v, q);
    C(1) = LognormalEndsConstC(v, q, a, 1);
    C(3) = LognormalEndsConstC(v, q, a, 3);
end
F = (C(1) + 1/a - lam1)^(3/2) - (C(1) + 1/a - lam_max)^(3/2);
F = F / sqrt(q^2 *C(3) + 1/a^3) * 2/3;
