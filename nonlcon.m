function [C, Ceq] = nonlcon(X)
minv = 8.0e-2;
maxv = 1.0e-1;
ep = 1.0e-3;

mmt = johnson_su_moments12(X(1:4));

% The variance of the innovations must be in [8.0e-2, 1.0e-1];
% The mean of the innovations must be within [-ep, ep]
C(1) = minv + ep - mmt(2);
C(2) = mmt(2) - maxv + ep;
C(3) = mmt(1) - ep;
C(4) = -ep - mmt(1);

Ceq = 0;
