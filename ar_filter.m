% Apply the operator
function y = ar_filter(x, C)
n = length(C) - 1;
l = length(x);
y = NaN(l - n, 1);

for m = l:-1:n+1
    y(m - n) = C * x(m:-1:m-n);
end
