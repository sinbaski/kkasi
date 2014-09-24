function [f, fder, hessian] = LognormalEndsFun2(X, v, q, s)
f = sign(s) * getfield(LognormalBlue(X(1) + i*X(2), v, q), {1});
fder = sign(s) * getfield(LognormalBlueDer(X(1) + i*X(2), v, q),...
                          {1; 1:2})';
hessian = NaN(2, 2);
A = LognormalBlueDer2(X(1)+i*X(2), v, q);
B = LognormalBlueDer2(X(1)-i*X(2), v, q);
hessian(1, 1) = (A+B)./2;
hessian(1, 2) = (A-B).*i./2;
hessian(2, 1) = hessian(1, 2);
hessian(2, 2) = -hessian(1, 1);
hessian = hessian .* sign(s);





