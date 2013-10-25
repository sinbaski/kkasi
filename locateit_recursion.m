% X: the vector
% pos: the starting position in X
% func: how to compare elements of X.
% l: expected distance
% n: recursion counter
% s: step size of the last recursion
function index = locateit_recursion(X, pos, l, func, param, n, s)
if n == 0
    s = [l, l/exp(1)] ./ exp(1);
    stepsize = int32(s(1));
elseif n == 1
    stepsize = int32(s(2));
else
    s = s./exp(1);
    stepsize = int32(s(2));
end
if stepsize == 0
    stepsize = 1;
end

odd = mod(n, 2);
if odd; stepsize = -stepsize; end

index = pos;
cond = 1;
while cond
    index = index + stepsize;
    if index > size(X, 1)
        index = size(X, 1);
    elseif index < 1
        index = 1;
    end
    %    fprintf('n=%d, pos=%d, index=%d\n', n, pos, index);
    cond = ((~odd && func(X, param(1), index, param(2)) < 0) ||...
            (odd && func(X, param(1), index, param(2)) > 0));
    if cond && index == size(X, 1)
        index = NaN;
        return;
    end
end

if func(X, param(1), index, param(2)) == 0
    return;
else
    index = locateit_recursion(X, index, l, func, param, n+1, s);
end
