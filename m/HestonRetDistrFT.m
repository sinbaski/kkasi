function ft = HestonRetDistrFT(model, p, t)
rho = model.Correlation;
rho = rho(2, 1);
mu = model.Return;
theta = model.Level;
gamma = model.Speed;
kappa = model.Volatility;

Gamma = gamma +i*rho*kappa*p;
Omega = sqrt(Gamma.^2 + kappa^2 * (p.^2 - i*p));
ft = gamma*theta/kappa^2 * Gamma*t ...
     -2*gamma*theta/kappa^2 * log(...
         cosh(t*Omega./2) + (Omega.^2 - Gamma.^2 + 2*gamma*Gamma)./ ...
         (2*gamma.*Omega).*sinh(Omega*t/2)...
         );

         
         