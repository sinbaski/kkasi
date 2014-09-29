% v: scalar
% q: scalar
% a: scalar or column vector
% return: scalar or column vector
function [f, fder] = LognormalEndsFun2(v, q, a)
[f, fder] = LognormalEndsFun1(v, q, a);
f = a - f^(-1/2);
fder = 1 + (1/2)*f.^(-3/2).*fder;
