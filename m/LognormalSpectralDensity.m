clear all
% close all
mysql = get_mysql();
name='SP500';
stmt = sprintf('select symbol from %s_components;', name);
symbols = fetch(mysql, stmt);

N = size(symbols, 1);
numRec = NaN(N, 1);
day1 = '2005-01-02';
day2 = '2015-01-23';
for k = 1 : N
    data = fetch(mysql, sprintf(['select count(*) ' ...
                        'from %s_US where day between "%s" and "%s";'], ...
                                strrep(symbols{k}, '.', '_'), day1, ...
                                day2));
    numRec(k) = cell2mat(data(1, 1));
end
T = max(numRec) - 1;
to_include = find(numRec == T + 1);
% d1 = datenum(day1);
% d2 = datenum(day2);
% days = ones(1, d2 - d1 + 1);
% for k = to_include'
%     stmt = sprintf('select distinct day from %s_US;', strrep(symbols{k}, '.', '_'));
%     data = datenum(fetch(mysql, stmt));
%     data >= d1 & data <= d2
% end

N = length(to_include);
R = NaN(T, N);
V = NaN(N, 1);
j = 1;
gam = -psi(1);

for k = to_include'
    % stmt = sprintf(['select closing from %s_US ' ...
    %                 'where day between "%s" and "%s" order by day;'], ...
    %                strrep(symbols{k}, '.', '_'), day1, day2);
    stmt = sprintf(['select closing, (high-low)/low from %s_US ' ...
                    'where day between "%s" and "%s" order by day;'], ...
                   strrep(symbols{k}, '.', '_'), day1, day2);
    data = cell2mat(fetch(mysql, stmt));
    R(:, j) = price2ret(data(:, 1));
    A = log(data(:, 2));
    V(j) = var(A);
    mu = mean(A);
    % a = kurtosis(R(:, j));
    % b = var(R(:, j));
    % V(j) = log(a/3)/4;
    % mu = log(b/sqrt(a/3))/2;
    R(:, j) = R(:, j) .* exp(-mu);
    j = j + 1;
end
close(mysql);

C = R' * R / T;
ev = sort(eig(C), 'descend');
ev = ev(11:end) ./ sum(ev(11:end)) .* sum(ev);

% n = floor(N/q);
% ev = NaN(N, 1);
% C = [];
% for offset = 0 : T - n
%     A = R(1:n, :);
%     B = R(1+offset : n+offset, :);
%     D = A' * B / n;
%     if isempty(C)
%         C = A./(T - n + 1);
%     else
%         C += A./(T - n + 1);
%     end
% end
% sections = floor(T/n);
% ev = NaN(N * sections, 1);
% for k = 1 : sections
%     A = R(n*(k-1)+1 : n*k, :);
%     C = A' * A ./ n;
%     ev(N*(k-1)+1 : N*k) = eig(C);
% end

% ev = sort(ev, 'descend');
% normalizer = ev(1);
% ev = ev ./ normalizer;

% if exist(sprintf('%s_logvol_variance_q%.2f.mat', name, q), 'file') == 2
%     load(sprintf('%s_logvol_variance_q%.2f.mat', name, q), 'vhat');
% else
%     ops = statset('mlecustom');
%     ops.MaxFunEvals = 600;
%     ops.MaxIter = 300;
%     vhat = mle(ev, 'pdf', @(x, va) abs(normalizer*imag(LognormalGreen(normalizer*x, va, q))/pi), ...
%                'start', mean(V), 'lowerbound', min(V), 'upperbound', max(V), ...
%                'options', ops);
%     save(sprintf('%s_logvol_variance_q%.2f.mat', name, q), 'vhat');    
% end

[X, Y] = epdf(ev, 1, min(ev), max(ev), 100, '');
% sample = ev;
% save('/tmp/DataSample.mat', 'sample');
% !Rscript ../r/EstimateDensity.r
% load('/tmp/DensityFunction.mat', 'density');
% X = density(:, 1);
% Y = density(:, 2);

% [X, Y] = epdf(ev, 1, min(ev), max(ev), 400, '');
% plot(X, Y, ev, density, 'LineWidth', 2);
% title(sprintf('q = %.2f', q));
% grid on
% vhat = mean(V);

q = N/T;
v = mean(V);
G = LognormalGreen(X, v, q);
density = -imag(G)/pi;

MPdensity = MarcenkoPasturPDF(X, [q, exp(v/2)]);

plot(X, Y, X, density, X, MPdensity, 'LineWidth', 2);
title(sprintf('S&P500, N=442, T=2531, q = %.2f, ', q));
grid on
legend('Empirical', 'Lognormal', 'MP');

% dx = X(2) - X(1);

% %% Hellinger
% % Hellinger distance between the empirical and the theoretical
% H1 = sum((sqrt(density) - sqrt(Y)).^2 * dx)^(1/2) / sqrt(2);
% % Hellinger distance between the empirical and the theoretical
% H2 = sum((sqrt(MPdensity) - sqrt(Y)).^2 * dx)^(1/2) / sqrt(2);

% %% Kullback-Leibler divergence of the theoretical from the empirical
% I = Y > 0;
% KL1 = sum(log(Y(I)./density(I)) .* Y(I) * dx);
% KL2 = sum(log(Y(I)./MPdensity(I)) .* Y(I) * dx);
