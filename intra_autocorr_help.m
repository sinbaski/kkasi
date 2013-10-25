clear all
load('intra_autocorr.mat', 'data');
cen = linspace(min(data.delta_t), max(data.delta_t), 20);
dt = cen(2) - cen(1);

X = ones(1, length(cen)) * NaN;
for n = 1:length(cen)
    I = data.delta_t >= cen(n) - dt/2 & ...
        data.delta_t < cen(n) + dt/2;
    X(n) = mean(data.auto_corr(I));
end
plot(cen, X);
