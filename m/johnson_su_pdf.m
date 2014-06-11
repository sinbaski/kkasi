%% johnson_su: a structure containing the parameters:
% Gamma, Delta, Xi, Lambda

function y = johnson_su_pdf(johnson_su, x)
%% The probability density function of the Johnson Su distr.
% x:  vector or scalar. The point(s) at which the pdf is evaluated.

Gamma = johnson_su.Gamma;
Delta = johnson_su.Delta;
Xi = johnson_su.Xi;
Lambda = johnson_su.Lambda;

u = (x - Xi) ./ Lambda;

a = 1/sqrt(2*pi);

b1 = asinh(u) .* Delta + Gamma;
b = exp(-b1.^2/2);

c1 = (Delta/Lambda);
c2 = (1 + u.^2).^(-1/2);
c = c1 .* c2;

y = a.*b.*c;
