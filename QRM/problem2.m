clear all
close all
L1 = importdata('DAX.txt', '\t');
L1 = -L1(:, 2);

%% compound into 10-day returns
L10 = NaN(length(L1)/10, 1);
for j = 1:length(L1)/10
    L10(j) = sum(L1((j-1)*10+1 : j*10));
end

L = L10;

n = length(L);
alfa = 0.99;
k = ceil(n*(1 - alfa));
Ls = sort(L, 'descend');
Ls(k);

fprintf('The VaR is the %dth upper order statistic.\n', k);

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
if m > 0
    m = m - 1;
end

fprintf('The upper bound is the %dth upper order statistic.\n', m+1);
