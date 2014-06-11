function mmt = johnson_su_moments(args)

mmt = NaN(1, 4);

mmt(1:2) = johnson_su_moments12(args);
mmt(3:4) = johnson_su_moments34(args(1:2));
