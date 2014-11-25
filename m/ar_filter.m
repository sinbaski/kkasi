% Apply the operator
function y = ar_filter(C, x)
n = length(C);
l = length(x);
y = NaN(l - n, 1);

for m = l:-1:n+1
    y(m - n) = x(m) - C * x(m-1:-1:m-n);
end
