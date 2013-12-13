function y = RealWishart_P(ain, x)
persistent a1 a n ya k u1 u2;

if isempty(a1) || ain ~= a1
    a1 = ain;
    ya = NaN(size(x));
    if mod(a1*2, 2) == 0
        n = a1;
        a = 0;
        ya = ones(size(x));
    else
        n = a1 - 1/2;
        a = 1/2;
        ya = 2*normcdf(sqrt(2*x)) - 1;
    end
    k = [0:n-1];
    u1 = factorial(k);
    u2 = NaN(size(k));
    for m = 1:length(k)
        u2(m) = sqrt(pi) * prod(1:2:(2*k(m)+1)) / 2^(k(m)+1);
    end
end

y1 = NaN(size(x));

if a == 0
    for j = 1 : length(x)
        y1(j) = sum(x(j).^(a+k) ./ u1);
    end
else
    for j = 1 : length(x)
        y1(j) = sum(x(j).^(a+k) ./ u2);
    end
end
y = ya - exp(-x) .* y1;
