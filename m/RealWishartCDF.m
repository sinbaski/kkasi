function p = RealWishartCDF(x, N, T)
% Gamma = @(m, a) pi^(m*(m-1)/4) .* prod(gamma(a - ([1:m] - 1)./2));
P = @RealWishart_P;

% I = @(a, b, x1) 1./gamma(a) .* integral(@(t) t.^(a-1) .* exp(-t) .* ...
%                                        gammainc(b, t), 0, x1);

nmax = max(N, T);
nmin = min(N, T);
if mod(nmin, 2) == 1
    nmat = nmin + 1;
else
    nmat = nmin;
end
alpha = (nmax - nmin - 1) / 2;

% K = pi^(nmin^2/2);
% K = K/2^(nmin*nmax/2);
% K = K/Gamma(nmin, nmax/2);
% K = K/Gamma(nmin, nmin/2);
% K1 = K*2^(alpha*nmat + nmat*(nmat+1)/2) * prod(alpha + [1:nmat]);
[nf, df] = RealWishart_K1(nmin, nmax);
A = zeros(nmat, nmat, length(x));
b = NaN(length(x), nmin);
for u = 1:nmin
    b(:, u) = (P(alpha+u, x/2)).^2 / 2;
    m = 2*(alpha-1/2) + 2*u;
    for v = u : nmin-1
        l = v - u;
        if mod(2*alpha, 2) == 1 && mod(l, 2) == 1
            k = 1:m/2;
            s1 = prod(2*k ./ (2*k-1));

            k = 1:(l+1)/2;
            s2 = prod((m+2*k-2) ./ (m+l+2*k));
            
            s = s1 * s2 /(m*pi);
        elseif mod(2*alpha, 2) == 1 && mod(l, 2) == 0
            k = 1:m/2;
            s1 = prod(2*k ./ (2*k-1));
            
            if l == 0
                s2 = 1;
            else
                k = 1 : l/2;
                s2 = prod((m+2*k) ./ (m+l+2*k+1));
            end
            s = s1 * s2 / (m+l+1) / pi;
        elseif mod(2*alpha, 2) == 0 && mod(l, 2) == 1
            k = 1:m/2-1;
            s1 = prod(2*k ./ (2*k-1));
            
            k = 1:(l+1)/2;
            s2 = prod((m+2*k-3) ./ (m+l+2*k-1));
            s = s1 * s2 / 2;
        elseif mod(2*alpha, 2) == 0 && mod(l, 2) == 0
            k = 1:m/2-1;
            s1 = prod(2*k ./ (2*k-1));
            
            k = 1:l/2+1;
            s2 = prod((m+2*k-3) ./ (m+l+2*k-2));
            s = s1 * s2 / 2;
        end
        b(:, v+1) = b(:, v) - s .* P(2*alpha + (u+v), x);
        A(u, v+1, :) = ...
            P(alpha+u, x/2) .* P(alpha+v+1, x/2)...
            - 2*b(:, v+1);
    end
end
if mod(nmin, 2) == 1
    for u = 1:nmin
        A(u, nmat, :) = P(alpha + u, x/2) ./ 2^(alpha + u) ./ ...
            gamma(alpha + nmin + 1);
    end
    for v = 1:nmin
        A(nmat, v, :) = -A(v, nmat, :);
    end
    A(nmat, nmat, :) = 0;
end

pf = NaN(1, length(x));
n = length(nf);
for k = 1:length(x)
    A(:, :, k) = A(:, :, k) - A(:, :, k)';
    A(:, :, k) = A(:, :, k) .* prod((nf./df).^(1/nmat));
    pf(k) = sqrt(det(A(:, :, k)));
    if pf(k) ~= 0
        p(k) = prod(nf .* pf(k)^(1/2/n) ./ df);
    else
        p(k) = 0;
        continue;
    end
end
