clear all
close all
L = importdata('DAX.txt', '\t');
L = -L(:, 2);
n = length(L);
alfa = 0.99;
k = ceil(n*(1 - alfa));
Ls = sort(L, 'descend');
Ls(k);

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

fprintf('The lower bound is the %dth upper order statistic.\n', m);

accum = alfa^n;
bita = 5.0e-2;
m = 0;
while accum <= bita/2
    m = m + 1;
    A = 1:m;
    B = n-m+1:n;
    accum = accum + prod(B ./ A .* (1 - alfa)) * alfa^(n-m);
end
m = m - 1;
