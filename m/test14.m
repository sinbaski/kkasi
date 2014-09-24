clear all
close all
sig = 0.5;
phi = 0;
q = 0.1;
v = sig^2 ./ (1 - phi^2);
a = erf(sqrt(2*v));
b = 1 - (1/2) * exp(2*v) * erf(sqrt(8*v));

load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/Eig-' ...
              'sig%.4f-phi%.4f.mat'],sig, phi), 'ev');

L = 300;
ev = reshape(ev, 1, prod(size(ev)));
M = max(ev);
lam = linspace(M*0.3, M, L+1);
% y1 = hist(ev(ev > b*0.3), lam);
% y1 = y1(1:end-1)./length(ev)./(lam(2)-lam(1));
lam = lam(1:end-1);

% y1 = y1./(lam(2)-lam(1));

G = LognormalGreen(lam, v, q);
x = real(G);
y = imag(G);

lam_com=exp(2*v) + exp(4*v)*sqrt(-q.*y) + sqrt(q)*exp(4*v)./sqrt(-y);
lam2 = sqrt(q)*exp(4*v)./sqrt(-y);

% lam_com=(exp(-8*v)*(8*exp(4*v)*v*(exp(4*v)-sqrt(q*y))+(exp(8*v)*(-1+8*v)-q*y+2*exp(4*v)*sqrt(q*y)).*log(exp(4*v)*(1+4*v)-sqrt(q*y))+(exp(8*v)*(1-8*v)+q*y-2*exp(4*v)*sqrt(q*y)).*log(exp(4*v)*(-1+4*v)+sqrt(q*y))))/(16*sqrt(2*pi)*v^(3/2))+sqrt(q*exp(8*v)./y);

% lam_com = real(lam_com);

% Y=q*exp(8*v)./(lam-exp(2*v)).^2;
% Y = -Y.*pi;
