function [F, FDer] = LognormalGreenFun1(Z, v, q, B)
F = LognormalBlue(Z, v, q)' - B;
FDer = LognormalBlueDer(Z, v, q);
return;
