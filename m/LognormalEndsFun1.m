function [c, ceq] = LognormalEndsFun1(X, v, q, lb)
ceq = getfield(LognormalBlue(X(1) + i*X(2), v, q), {2});
c = lb - getfield(LognormalBlue(X(1) + i*X(2), v, q), {1});

