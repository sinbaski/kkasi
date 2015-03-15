clear all
close all

N = 600;
fold = 1.0e+4;
q = 1.0;
sig = sqrt(1/2);
phi = 0;

load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eigvec-' ...
              'sig%.4f-phi%.4f.mat'], q, sig, phi), 'eigvec');
npr = @(x) 1/sum(x.^4)/N;
NPR = NaN(2, fold);
for k = 1 : fold
    NPR(1, k) = npr(eigvec(:, k, 1));
    NPR(2, k) = npr(eigvec(:, k, 2));
end
sample = NPR(1, :);
save('/tmp/DataSample.mat', 'sample');
!Rscript ../r/EstimateDensity.r
load('/tmp/DensityFunction.mat', 'density');
X = density(:, 1);
Y = density(:, 2);
plot(X, Y, 'b-');
hold on
sample = NPR(2, :);
save('/tmp/DataSample.mat', 'sample');
!Rscript ../r/EstimateDensity.r
load('/tmp/DensityFunction.mat', 'density');
X = density(:, 1);
Y = density(:, 2);
plot(X, Y, 'g-');

title('PDF of the Normalized Participation Ratio');
legend('X_1', 'X_{300}');
grid on




% load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
%               'sig%.4f-phi%.4f.mat'], q, sig, phi), 'ev');
% load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Diagelem-' ...
%               'sig%.4f-phi%.4f.mat'], q, sig, phi), 'diagelem');
% [yev, xev] = ecdf(ev(1, :));
% [yelem, xelem] = ecdf(diagelem(1, :));
% Qev = quantile(ev(1, :), [0.98, 0.999]);
% Qelem = quantile(diagelem(1, :), [0.98, 0.999]);
% loglog(xev, 1 - yev, xelem, 1 - yelem);
% % xlim([min([Qev, Qelem]), max([Qev, Qelem])]);
% ylim([1.0e-3, 1]);
% legend('\lambda_{max}', 'C_{ii}');
% grid on
% title('Tail function of C_{ii} and \lambda_{max}');

% X = sort(diagelem(1, :), 'descend');
% tailIndex = NaN(1, length(X)-1);
% for k = 2 : length(X)
%     tailIndex(k-1) = 1/mean(log(X(1:k-1) ./ X(k)));
% end
% plot(1:length(X)-1, tailIndex);
% xlabel('#Upper order statistics');
% ylabel('Tail index estimate');
% xlim([1, 400]);
% ylim([2, 5]);
% grid on
