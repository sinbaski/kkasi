clear all
close all
dof = 3;
d = (4-dof)/(2*(dof - 1)) * 0.9;
p = 10;
n = floor(p^(1/d));

specs1 = {'b', 'g', 'c', 'r'};

load(sprintf(['~/matfys/data/student3/Eig.mat']), 'ev');
fold = size(ev, 2);
ev = sort(ev, 'descend');

E = exprnd(1, p, fold);
Gam = cumsum(E, 1);

% for a = 1:4
%     X = ev(a+1, :)./ev(a, :);
%     Y = (Gam(a+1, :) ./ Gam(a, :)).^(-2/dof);
    
%     subplot(2, 2, a);
%     qqplot(X, Y);
%     ylim([0, 2]);
%     grid on
%     title(sprintf('\\alpha = %.1f, p=%d, n=%d', dof, p, n));
%     xlabel(sprintf('\\lambda_{(%d)}/\\lambda_{(%d)}', a+1, a));
%     ylabel(sprintf('(\\Gamma_{%d}/\\Gamma_{%d})^{-2/\\alpha}', a+1, a));
% end
hold on
for a = 1:4
    egenvarde = ev(a+1, :)./ev(a, :);
    
    % subplot(2, 2, a);
    [F, x] = ecdf(egenvarde);
    plot(log10(x), log10(F), specs1{a});
    grid on    
    title(sprintf('\\alpha = %.1f, p=%d, n=%d', dof, p, n));
    ylabel('lg(F)');
    texts{a} = sprintf('\\lambda_{(%d)}/\\lambda_{(%d)}', a+1, a);
end
hold off
legend(texts);
ylim([-3, 0]);
xlabel(sprintf('$\\lg(\\lambda_{(i+1)}/\\lambda_{(i)})$', a+1, ...
               a), 'Interpreter', 'Latex', 'Fontsize', 18);


