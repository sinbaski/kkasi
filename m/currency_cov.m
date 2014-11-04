clear all
close all
mysql = get_mysql();
p = 5;
n = 250;
periods = 36;
currencies = {'GBP', 'CHF', 'SEK', 'DKK', 'JPY'};
data = fetch(mysql, sprintf(['select gbp,chf,sek,dkk,jpy from USD_Rates where ',...
                    'gbp is not null and ', ...
                    'chf is not null and ', ...
                    'sek is not null and ', ...
                    'dkk is not null and ', ...
                    'jpy is not null order by day desc limit %d;', ...
                    ], p*n*periods+1));
% data = fetch(mysql, sprintf(['select gbp, jpy from USD_Rates where ',...
%                     'gbp is not null and ', ...
%                     'jpy is not null order by day desc limit %d;', ...
%                     ], p*n*periods+1));
close(mysql);
rates = 1./flipud(cell2mat(data));
X = price2ret(rates);
EC = cov(X);

lambda = NaN(5, periods);
for c = 0: periods-1
    lambda(:, c+1) = sort(eig(cov(X(c*n+1 : (c+1)*n, :)) - EC), 'descend');
end
I = lambda(1, :) > 0 & lambda(2, :) > 0;
A = log10(lambda(1, I)./lambda(2, I));
A = A./std(A);

B = log10(1./rand(1, 1000));
B = B./std(B);
qqplot(A, B);

xlabel(['Quantile of $X/std(X)$, where X = $\log(\lambda_{(1)}/\' ...
        'lambda_{(2)})$'], 'Interpreter', 'Latex', 'Fontsize', 14);
ylabel('Quantile of $Y/std(Y)$, where Y = $log(1/U)$', 'Interpreter', ...
       'Latex', 'Fontsize', 14);
grid on
title(['QQ-plot In cases \lambda_{(1)} > \lambda_{(2)} > 0.   '...
      'p = 5, n = 250, 36 periods']);

figure;
for c = 1 : size(X, 2)
    subplot(2, 3, c);
    autocorr(X(:, c), 10);
    title(sprintf('ACF of %s prices in USD', currencies{c}));
    ylim([-0.05, 0.2]);
end

xcorr = corr(X);
fd = fopen('xcorr.txt', 'wt');
fprintf(fd, '% 6s', ' ');
for c = 1:p
    fprintf(fd, '% 6s', currencies{c});
end
fprintf(fd, '\n');
for c = 1:p
    fprintf(fd, '% 6s', currencies{c});
    for d = 1 : p
        fprintf(fd, '% 6.2f', xcorr(c,d));
    end
    fprintf(fd, '\n');
end
fclose(fd);

% gbp = X(:, 1);
% tail = gbp(gbp > 0.01);
% [y, x] = ecdf(tail);
% y = 1 - y;
% I = y > 10^(-1.2);
% x = x(I);
% y = y(I);
% plot(log10(x), log10(y), 'bx');
% grid on

% P = polyfit(log10(x), log10(y), 1);
% alfa = -P(1);
% d = (4-alfa)/(2*(alfa - 1));
% p^(1/d)






