% spectral density of the Wishart-Levy matrix
clear all
generate_data = 0;

n = 500; % number of time series
T = 2000;
%distr.name = 'lorentzian';
distr.name = 'student t';
distr.prmt = 3;

tau = 1;

%auto_corr = @(t) exp(-t/tau);
auto_corr = @(t) 0;

lb = NaN;
ub = NaN;
total_num = 0;

if generate_data
    for l = 1:2000
        filename = sprintf('data/StudentWishart/eig_%d.mat', l);
        
        R = gen_ret_mtx(n, T, distr, auto_corr);
        %C = cpt_cov(R, distr);
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
        
        save(filename, 'lambda');
        fprintf('Saved file %s\n', filename);
    end
    save('data/StudentWishart/lbub.mat', 'lb', 'ub', 'total_num');
end

% load('data/StudentWishart/lbub.mat');
% cen = linspace(lb, ub, 1000);
% dx = cen(2) - cen(1);
% stats = zeros(1, length(cen));

eig = [];
for l = 1:2000
    filename = sprintf('data/StudentWishart/eig_%d.mat', l);
    load(filename, 'lambda');
    eig = [eig; lambda(lambda < 20)];
    %plot(lambda, ones(length(lambda), 1), 'o');
    %stats = stats + hist(lambda, cen);
end
cen = linspace(min(eig), max(eig), 1000);
dx = cen(2) - cen(1);
empirical = double(hist(eig, cen)) ./ (length(eig) * dx);
%empirical = double(stats) ./ (total_num * dx);
plot(cen, empirical, 'bo');
