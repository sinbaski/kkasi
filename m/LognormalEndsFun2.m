% v: scalar
% q: scalar
% a: scalar or column vector
% return: scalar or column vector
% function [f, fder] = LognormalEndsFun2(v, q, u)
function f = LognormalEndsFun2(v, q, u)
% [f, fder] = LognormalEndsFun1(v, q, log(u));
f = LognormalEndsFun1(v, q, log(u));
f = log(u)^(-2) - f;
% fder = (-2*log(u)^(-3) - fder)/u;
[log(u), f]
;
