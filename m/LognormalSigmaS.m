% Moment generating function of the matrix
% /                 \
% |  s_1            |
% |       s_2       |
% |          ,      |
% \            s_T  /
% Z1: column vector, complex values
% v1: real positive
% ret: column vector, complex values
function ret = LognormalSigmaS(Z1, v1)
persistent S Z v;

if (isempty(Z) || length(Z) ~= length(Z1) ||...
    sum(Z == Z1) ~= length(Z) || ...
    isempty(v) || v1 ~= v)
    v = v1;
    Z = Z1;
else
    ret = S;
    return;
end

for n = 1 : length(Z)
    lsqnonlin(@(X) LognormalSigmaM(X(1) + i*X(2), v1) - Z(n)
end
%% First find the N transform as the invserse of M
