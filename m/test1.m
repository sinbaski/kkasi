clear all
close all

mysql = get_mysql();
stmt = sprintf('select symbol from DAX_components;');
data = fetch(mysql, stmt);
p = 29;
T = 3000;

R = NaN(p, T);
counts = NaN(1, p);
for k = 1 : p
    closing = fetch(mysql, sprintf(['select closing from %s order by day desc ' ...
                        'limit %d;'], strrep(data{k}, '.', '_'), ...
          T+1));
    R(k, :) = price2ret(cell2mat((flipud(closing))))';
end
close(mysql);

EC = cov(R');

NbrPeriods = 10;
sampleSize = T/NbrPeriods;
ev = NaN(p, NbrPeriods);

for m = 1 : NbrPeriods
    A = R(:, (m-1)*sampleSize+1 : m*sampleSize);
    ev(:, m) = reshape(sort(eig(A*A' - EC), 'descend'), p, 1);
    % ev(:, m) = reshape(sort(eig(A*A'), 'descend'), p, 1);
end

ratios = NaN(5, NbrPeriods);
for k = 1:5
    ratios(k, :) = k*log((ev(k, :)./ev(k+1, :)));
end
A = reshape(ratios, 1, numel(ratios));

% qqplot(A./std(A), ProbDistUnivParam('exponential', 1));
% ylabel('Quantiles of $k[\ln{\lambda_{(k)} - \ln\lambda_{(k+1)}}]$', ...
%        'Interpreter', 'Latex', 'Fontsize', 14);
% grid on

[y, x] = ecdf(A);
plot(x(1:end-1), log(1 - y(1:end-1)), '+');
P = polyfit(x(1:end-1), log(1 - y(1:end-1)), 1);
grid on
title(['Tail function of $k \ln{\lambda_{(k)} \over \lambda_{(k+1)}}$ ' ...
       '$k=1,\dots,5$ on semi-log scale'], 'Interpreter', 'Latex', 'Fontsize', 16);
ylabel('$\ln\bar F(x)$', 'Interpreter', 'Latex', 'Fontsize', 14);
xlabel('x');

% results = NaN(1, p-1);
% U = sort(rand(500, p-1), 'descend');
% for k = 1: p-1
%     eig_sample = log(lambda(p, :)./lambda(p-k, :));
%     eig_sample = eig_sample./std(eig_sample);
%     uni_sample = log(U(k, :));
%     uni_sample = uni_sample ./ std(uni_sample);
%     results(k) = kstest2(eig_sample, uni_sample);
% end

% fd = fopen('results.txt', 'wt');
% for k = 1 : p - 1
%     fprintf(fd, '%d, ', results(k));
%     fprintf(fd, '\n');
% end
% fclose(fd);

% Check whether lambda()
% C = R*R';
% lambda = flipud(sort(eig(C)));
% lambda = lambda ./ lambda(1);


