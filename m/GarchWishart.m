% spectral density of the Wishart-Levy matrix
clear all
generate_data = 0;

n = 500; % number of time series
T = 2000;
%distr.name = 'lorentzian';
distr.name = 'garch';
distr.prmt = struct('Name', 'Gaussian');

tau = 1;

%auto_corr = @(t) exp(-t/tau);
auto_corr = @(t) 0;

lb = NaN;
ub = NaN;
total_num = 0;

if generate_data
    for l = 1:2000
        eig_file = sprintf('data/GarchWishart/eig_%d.mat', l);
        ret_file = sprintf('data/GarchWishart/ret_%d.mat', l);
        
        R = gen_ret_mtx(n, T, distr, auto_corr);
        save(ret_file, 'R');

        sigma = std(R');
        R = diag(1./sigma) * R;
        
        C = (1/T) * R * R';
        lambda = eig(C);
        
        a = min(lambda);
        b = max(lambda);
        total_num = total_num + length(lambda);
        
        if isnan(lb) || lb > min(lambda)
            lb = a;
        end
        if isnan(ub) || ub < max(lambda)
            ub = b;
        end
        
        save(eig_file, 'lambda');
        fprintf('Saved file %s\n', eig_file);
    end
    save('data/GarchWishart/lbub.mat', 'lb', 'ub', 'total_num');
end

% load('data/GarchWishart/lbub.mat');
% cen = linspace(lb, ub, 400);
% dx = cen(2) - cen(1);
% stats = zeros(1, length(cen));

eig = [];
for l = 1:2000
    eig_file = sprintf('data/GarchWishart/eig_%d.mat', l);
    load(eig_file, 'lambda');
    eig = [eig; lambda(lambda < 10)];
    %plot(lambda, ones(length(lambda), 1), 'o');
    %stats = stats + hist(lambda, cen);
end

b = max(eig) * 1.1;
a = min(eig) * 0.9;
q0 = (sqrt(b) - sqrt(a))^2/(sqrt(a) + sqrt(b))^2;
var0 = (a + b) / (2 + 2*q0);

cen = linspace(a, b, 500);
dx = cen(2) - cen(1);
empirical = double(hist(eig, cen)) ./ (length(eig) * dx);
[prmt, resnorm] = lsqcurvefit(@MarcenkoPasturPDF, [q0, var0], cen, empirical);
theoretical = MarcenkoPasturPDF(prmt, cen);
plot(cen, empirical, 'bo', cen, theoretical, 'r-');

