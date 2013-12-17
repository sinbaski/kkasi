function y = RealWishart_P(ain, x)
persistent a1 a n ya k u1 u2;

% if isempty(a1) || ain ~= a1
%     a1 = ain;
%     ya = NaN(size(x));
%     if mod(a1*2, 2) == 0
%         n = a1;
%         a = 0;
%         ya = ones(size(x));
%     else
%         n = a1 - 1/2;
%         a = 1/2;
        
%     end
%     % k = [0:n-1];
%     % u1 = factorial(k);
%     % u2 = NaN(size(k));
%     % for m = 1:length(k)
%     %     u2(m) = sqrt(pi) * prod(1:2:(2*k(m)+1)) / 2^(k(m)+1);
%     % end
% end
y = NaN(size(x));
if a == 0
    for j = 1 : length(x)
        inc = 0;
        for k = 0:n-1
            inc = inc + prod(x(j)./[1:k]);
        end
        y(j) = 1 - inc / exp(x(j));
    end
else
    for j = 1 : length(x)
        inc = 0;
        for k = 0:n-1
            inc = inc + prod(x(j)./[1/2:1:k+1/2]);
        end
        y(j) = erf(sqrt(x(j))) - inc/exp(x(j))/sqrt(pi*x(j));
    end
end

