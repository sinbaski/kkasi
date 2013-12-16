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
    % u1 = factorial(k);
    % u2 = NaN(size(k));
    % for m = 1:length(k)
    %     u2(m) = sqrt(pi) * prod(1:2:(2*k(m)+1)) / 2^(k(m)+1);
    % end
end

y1 = NaN(size(x));

if a == 0
    for j = 1 : length(x)
        y1(j) = 0;
        l = 0;
        done = 0;
        while ~done
            inc = 0;
            for k = 0 : min(n-1, l)
                inc = inc + prod(x(j) ./ [1:k, 1:l-k]) * (-1)^(l-k);
            end
            y1(j) = y1(j) + inc;
            if abs(inc) > 0.05
                l = l + 1;
            else
                done = 1;
            end
        end
    end
else
    for j = 1 : length(x)
        y1(j) = 0;
        l = 0;
        done = 0;
        while ~done
            inc = 0;
            for k = 0:l
                inc = inc + prod(x(j) ./ [1:l, [1:2*k+1]./2]) * (-1)^l;
            end
            y1(j) = y1(j) + inc;
            if abs(inc) > 0.05
                l = l + 1;
            else
                done = 1;
            end
        end
        y1(j) = y1(j) /sqrt(x(j)) / sqrt(pi);
    end
end
y = ya - y1;
