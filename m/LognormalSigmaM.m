% Moment generating function of the matrix diag(s_1, s_2, ..., s_T)
% Z1: column vector, complex values
% v1: real positive
% ret: column vector, complex values
function ret = LognormalSigmaM(Z1, v1)
persistent M Z v;

if (isempty(Z) || length(Z) ~= length(Z1) ||...
    sum(Z == Z1) ~= length(Z) || ...
    isempty(v) || v1 ~= v)
    v = v1;
    Z = Z1;
else
    ret = M;
    return;
end

l = 2000;
s = linspace(-6*sqrt(v), 6*sqrt(v), l);

M1 = 1./(repmat(Z, 1, l) - repmat(exp(2*s), length(Z), 1))...
     ./sqrt(2*pi*v);
M2 = exp(-s'.^2./2) .* (s(2) - s(1));
M = Z.*(M1 * M2)-1;
ret = M;


