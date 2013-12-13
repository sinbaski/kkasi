clear all
close all
load('/tmp/models.mat');

T = 1e6;
[LV, E, V] = simulate(model0, T, 'numPaths', 1);
ret = exp(LV) .* randn(size(LV));
r = (ret - mean(ret))/std(ret);
[p, x] = ecdf(r);

%% Fit to Gaussian + power-law tails
% param = mle(r, 'pdf',...
%             @(x, xi1, xi2)...
%             NormalPowerLaw_pdf([xi1, xi2], x),...
%             'start', sqrt([1.4, 1.42]), ...
%             'lowerbound', [1.0, 1.0]);
% powerlaw = NormalPowerLaw_cdf(param, x);


%% Fit to Johnson Su dist.
% js = johnson_su_params([mean(r), var(r), skewness(r), ...
%                     kurtosis(r)]);
% johnson = johnson_su_cdf(johnson_su_struct(js), x);
% plot(x, pt, 'g');


% [p1, p2] = extract_xpnt(r, 0.046/std(ret));
[p1, p2] = extract_xpnt(r, 0.045/std(ret));
fprintf('Left tail exponent: %.4f\nRight tail exponent: %.4f\n', ...
        p1, p2);
