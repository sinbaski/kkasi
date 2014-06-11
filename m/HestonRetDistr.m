% Stationary distribution of the returns in the Hestion model
function prob_den = HestonRetDistr(model, t, x)
rho = model.Correlation(2, 1);

prob_den = zeros(length(x), 1);
for n = 1:length(x)
    func = @(p) exp(i*p*x(n) + HestonRetDistrFT(model, p, t));
    lim = 50;
    Q = integral(func, -lim, lim);
    diff = NaN;
    errtol = 1.0e-4;
    increment = 10;
    while isnan(diff) || abs(diff) > errtol 
        diff = integral(func, -lim - increment, -lim) +...
               integral(func, lim, lim + increment);
        lim = lim + increment;
        Q = Q + diff;
    end
    prob_den(n) = Q;
end
