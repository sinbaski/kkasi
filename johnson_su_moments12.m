function mmt = johnson_su_moments12(args)

Gamma = args(1);
Delta = args(2);
Lambda = args(3);
Xi = args(4);


w = exp(1 ./ Delta.^2);
Omega = Gamma ./ Delta;

m1 = -w.^(1/2) .* sinh(Omega) .* Lambda + Xi;
m2 = (1/2) * (w-1) .* (w.*cosh(2*Omega) + 1) .* Lambda^2;
mmt = [m1, m2];
