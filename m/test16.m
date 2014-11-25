clear all
close all
sig = 0.2;
q = 0.1;
v = sig^2;

a = LognormalEndsConstA(v, q);
C1 = LognormalEndsConstC(v, q, a, 1);
C3 = LognormalEndsConstC(v, q, a, 3);
[u, lam2] = LognormalEnds(v, q);

blk = 3*sqrt(C3*q^2 + 1/a^3)/2 + (C1 + 1/a - lam2)^(3/2);
lam1 = C1 + 1/a - blk^(2/3);
