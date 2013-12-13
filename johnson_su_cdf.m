function P = johnson_su_cdf(johnson, x)
% Evaluate Prob(u < x)
P = NaN(1, length(x));
for n = 1:length(x)
    P(n) = integral(@(u) johnson_su_pdf(johnson, u), -Inf, x(n));
end
