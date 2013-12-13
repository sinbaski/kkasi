function param = johnson_su_params(mmt)
param = NaN(1, 4);
[Gamma, Delta] = johnson_su_tables(mmt(3), mmt(4));
if isnan(Gamma) || isnan(Delta)
    throw(MException('johnson_su:johnson_su_params', 'Gamma or Delta is NaN'));
end
param([1:2]) = fsolve(@(x) johnson_su_moments34(x) - mmt(3:4), [Gamma, Delta]);
m = johnson_su_moments12([param([1:2]), 1, 0]);
param(3) = sqrt(mmt(2)/m(2));
param(4) = mmt(1) - m(1)*param(3);

