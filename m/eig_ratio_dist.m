clear all
close all
p = 500;
T = 500;

load(sprintf(['../matfys/data/student3/Eig.mat']), 'ev');
fold = size(ev, 2);
ev = sort(ev, 'descend');
X = NaN(1, fold);
Y = NaN(1, fold);
alfa = 3;
N = 50000;
% specs1 = {'b', 'g', 'c', 'r'};
% specs2 = {'b--', 'g--', 'c--', 'r--'};
if alfa > 2 && alfa < 4
    K = NaN(1, N);
    K(1) = gamma(1- 2/alfa);
    for k = 2:N
        K(k) = K(k-1) * (k-1-2/alfa)/(k-1);
    end
end

for n = 1
    for k = 1 : fold
        X(k) = ev(1, k) / sum(ev(:, k));
        E = exprnd(1, 1, N);
        Gam = cumsum(E);
        if alfa > 2
            Y(k) = Gam(1)^(-2/alfa) / sum(Gam.^(-2/alfa) - K, 2);
        else
            Y(k) = Gam(1)^(-2/alfa) / sum(Gam.^(-2/alfa), 2);
        end
    end
    handle = figure;
    qqplot(X, Y);
    xlim([0.01, 0.2]);
    hgexport(handle, 'student3_qq.pdf', hgexport('factorystyle'), 'Format', 'pdf');
end

% [F1, val1] = ecdf(X);
% [F2, val2] = ecdf(Y);
% handle = figure;
% plot(val1, log(F1), val2, log(F2));

% ylabel('CDF');
% xlabel('\lambda_{(1)}/trace');



