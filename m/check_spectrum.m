clear all

filename = 'data/wishart_eig.txt';
lambda = load(filename);

cen = linspace(min(lambda), max(lambda), 400);
dx = cen(2) - cen(1);
empirical = double(hist(lambda, cen)) ./ (length(lambda) * dx);

% [q, sigma]
%prmt0 = [0.25, 1];
b = max(lambda) * 1.1;
a = min(lambda) * 0.9;


q0 = (sqrt(b) - sqrt(a))^2/(sqrt(a) + sqrt(b))^2;
var0 = (a + b) / (2 + 2*q0);

%phat = mle(lambda, 'pdf', @MarcenkoPasturPDF, 'start', [q0,
%var0]);
[prmt, resnorm] = lsqcurvefit(@MarcenkoPasturPDF, [q0, var0], cen, empirical);

theoretical = MarcenkoPasturPDF(prmt, cen);
plot(cen, empirical, 'bo', cen, theoretical, 'r-');
