clear all

n = 500;
T = 2000;

% figure
% colors = {'b', 'g', 'black', 'm', 'r'};

filename = sprintf('data/GarchWishart/ret_%d.mat', 1);
load(filename, 'R');

R = reshape(R', 1, n*T);

cen = linspace(min(R), max(R), 1000);
dr = cen(2) - cen(1);
rho = double(hist(R, cen)) ./ (n*T * dr);

pd = fitdist(R', 'Normal');
rho_normalfit = pdf(pd, cen);

studentPdf = @(n, X)...
    gamma((n+1)/2) ./ (...
    sqrt(n*pi) * gamma(n/2) * (1 + X.^2 ./ n).^((n+1)/2));

[prmt, resnorm] = lsqcurvefit(studentPdf, [3], cen, rho);


%nu = mle(R', 'logpdf', studentLogPdf, 'start', 3);


plot(cen, rho, 'b+', cen, rho_normalfit, 'g', cen, studentPdf(prmt, ...
                                                  cen), 'r');

