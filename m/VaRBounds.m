function [lb, ub] = VaRBounds(alfa, bita, n)
accum = 0;
bita = 5.0e-2;
m = n;
while accum <= bita/2
    A = 1:m;
    B = n-m+1:n;
    accum = accum + prod(B ./ A .* (1 - alfa)) * alfa^(n-m);
    m = m - 1;
end
m = m + 2;

lb = m;

accum = alfa^n;
bita = 5.0e-2;
m = 0;
while accum <= bita/2
    m = m + 1;
    A = 1:m;
    B = n-m+1:n;
    accum = accum + prod(B ./ A .* (1 - alfa)) * alfa^(n-m);
end
if m > 0
    m = m - 1;
end

ub = m+1;

