% v: scalar
% q: scalar
% return: scalar or column vector
function [end1, end2] = LognormalEnds(v, q)
options = optimset('Jacobian','on');
% [xmin, fval1] = fsolve(@(a) LognormalEndsFun2(v, q, a), -1, options);
% [xmax, fval2] = fsolve(@(a) LognormalEndsFun2(v, q, a), 1, options);

% [end1, fval1] = fsolve(@(a) LognormalEndsFun2(v, q*exp(-2*v), a), 4, options);
% end1 = end1.^(-2)/q;
end1 = NaN;
[end2, fval2] = fsolve(@(a) LognormalEndsFun2(v, q*exp(-4*v), a), 0.1, options);
end2 = 1/end2^2/q/exp(2*v);




