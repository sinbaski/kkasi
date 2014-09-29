function a = LognormalEndsConstA(v, q)
options = optimset('Jacobian','on');
[a, fval2] = fsolve(@(x) LognormalEndsFun2(v, q, x), 1, options);
