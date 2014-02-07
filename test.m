clear all
s = 1;
p = 0.5;

x = linspace(-10, 10, 1000);
a = linspace(-10, 10, 1000);
[x, a] = meshgrid(x, a);

% z = exp(-a.^2./2./(1-p^2) + a.*p.*exp(-a.*s).*x./(1 - p^2));
z = exp(-(a.*p - exp(-a.*s).*x).^2 ./ 2 ./ (1 - p^2));
% surf(x, a, z);

% camlight('headlight');
contour(x, a, z);
colormap(cool);
xlabel('x');
ylabel('a');
