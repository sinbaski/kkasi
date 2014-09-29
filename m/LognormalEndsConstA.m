function a = LognormalEndsConstA(v, q, a0)
% options = optimset('Jacobian','on');
% [u, fval2] = fsolve(@(x) LognormalEndsFun2(v, q, x), exp(a0), ...
%                    options);
% [u, fval2] = fsolve(@(x) LognormalEndsFun2(v, q, x), exp(a0));
[u, fval2] = fsolve(@(x) LognormalEndsFun2(v, q, x), exp(a0));
a = log(u);
;


% if a0 > 1
%     options = optimset('Jacobian','on');
%     [u, fval2] = fsolve(@(x) LognormalEndsFun2(v, q, x), exp(a0), ...
%                         options);
%     a = log(u);
% else
%     [a, fval2] = fsolve(@(x) exp(8*v) + 2*q*exp(18*v)*x + 3*q^2* ...
%                         exp(32*v)*x^2 - 1/x^2, a0);
% end
