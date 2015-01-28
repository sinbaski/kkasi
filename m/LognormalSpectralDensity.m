clear all
close all
mysql = get_mysql();
name='SP500';
stmt = sprintf('select symbol from %s_components;', name);
symbols = fetch(mysql, stmt);

N = size(symbols, 1);
numRec = NaN(N, 1);
day1 = '2000-01-02';
day2 = '2015-01-23';
for k = 1 : N
    data = fetch(mysql, sprintf(['select count(*) ' ...
                        'from %s_US where day between "%s" and "%s";'], ...
                                strrep(symbols{k}, '.', '_'), day1, ...
                                day2));
    numRec(k) = cell2mat(data(1, 1));
end
T = 3788;
to_include = find(numRec == T);
% d1 = datenum(day1);
% d2 = datenum(day2);
% days = ones(1, d2 - d1 + 1);
% for k = to_include'
%     stmt = sprintf('select distinct day from %s_US;', strrep(symbols{k}, '.', '_'));
%     data = datenum(fetch(mysql, stmt));
%     data >= d1 & data <= d2
% end

N = length(to_include);
R = NaN(T-1, N);
j = 1;
counter = NaN(1, N);
for k = to_include'
    stmt = sprintf('select closing from %s_US where day between "%s" and "%s" order by day;', ...
                                   strrep(symbols{k}, '.', '_'), day1, day2);
    closing = fetch(mysql, stmt);
    R(:, j) = price2ret(cell2mat(closing));
    counter(j) = length(closing);
    j = j + 1;
end
close(mysql);

c = 1;
q = 0.8;

n = floor(N/q);
sections = floor(T/n);
ev = NaN(N * sections, 1);
for k = 1 : sections
    A = R(n*(k-1)+1 : n*k, :);
    C = A' * A;
    ev(N*(k-1)+1 : N*k) = eig(C);
end
ev = sort(ev, 'descend');
% Remove the largest from each section
ev = ev(sections + 1 : end);
sig = ev(1);
ev = ev ./ sig;
v = 0.8;

ops = statset('mlecustom');
ops.MaxFunEvals = 600;
ops.MaxIter = 300;
vhat = mle(ev, 'pdf', @(x, va) abs(sig*imag(LognormalGreen(sig*x, va, q))/pi), ...
           'start', v, 'lowerbound', 0.4, 'upperbound', 0.95, ...
           'options', ops);

[X, Y] = epdf(ev, 1, min(ev), max(ev), 100, '');
plot(X, Y, ev, density, 'LineWidth', 2);
title(sprintf('q = %.2f', q));
grid on

G = LognormalGreen(sig.*X, vhat, q);
density = -sig*imag(G)/pi;

MPdensity = MarcenkoPasturPDF(X, [q, v]);

dx = X(2) - X(1);
KL_entropy1 = sum(Y .* dx .* log(Y./density));
KL_entropy2 = sum(Y .* dx .* log(Y./MPdensity));

%% compute the Kullback-Leibler distance and the Hellinger distance

% [X, edos] = epdf(ev, 2, min(ev), max(ev), 60, 'b-');
% plot(X, edos, lam, density, 'LineWidth', 2);
% grid on

