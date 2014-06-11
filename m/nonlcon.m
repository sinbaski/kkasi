function [C, Ceq] = nonlcon(param, model, MALags, obs)
minv = 8.0e-2;
maxv = 1.0e-1;
ep = 1.0e-3;
s = 33;

c = 1;
for k = MALags
    model.MA{k} = param(c);
    c = c + 1;
end
y = infer(model, obs);
mmt = [mean(y), var(y), skewness(y), kurtosis(y)];

% The variance of the innovations must be in [8.0e-2, 1.0e-1];
% The mean of the innovations must be within [-ep, ep]
C(1) = minv + ep - mmt(2);
C(2) = mmt(2) - maxv + ep;
C(3) = mmt(1) - ep;
C(4) = -ep - mmt(1);
C(5) = mmt(3) - 0.5;
C(6) = -mmt(3);
C(7) = 3.6 - mmt(4);
C(8) = mmt(4) - 8.0;

Ceq = 0;
