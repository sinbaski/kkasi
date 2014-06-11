function Y = frct_dfrc(X, C)
if isempty(C)
    Y = NaN;
end
n = length(C) - 1;
l = length(X);
Y = NaN(l-n, 1);

for m=l:-1:n+1
    Y(m-n) = C * X(m:-1:m-n);
end
