% v: scalar
% q: scalar
% a: scalar or column vector
% return: scalar or column vector
function [f, fder] = LognormalEndsFun2(v, q, a)
[f, fder] = LognormalEndsFun1(v, q, a);
f = 1./a.^2 - f;
fder = -2./a.^3 - fder;
