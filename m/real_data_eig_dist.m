clear all
close all
mysql = get_mysql();
name='SP500';
stmt = sprintf('select symbol from %s_components;', name);
symbols = fetch(mysql, stmt);

p = size(symbols, 1);
T = 3000;
numRec = NaN(p, 1);

for k = 1 : p
    numRec(k) = cell2mat(fetch(mysql, sprintf(['select count(*) from %s_US'], ...
                                              strrep(symbols{k}, '.', '_'))));
end

to_include = find(numRec >= T+1);
p = length(to_include);

R = NaN(T, p);
c = 1;
for k = to_include'
    closing = fetch(mysql, sprintf(['select closing from %s_US order by day desc ' ...
                        'limit %d;'], strrep(symbols{k}, '.', '_'), ...
          T+1));
    R(:, c) = price2ret(cell2mat(flipud(closing)));
    R(:, c) = R(:, c) - mean(R(:, c));
    c = c + 1;
end
close(mysql);

k0 = 5;
M = NaN(p, p, k0+1);
for k = 1 : k0 + 1
    S = R(1:T-(k-1), :)' * R(k:T, :);
    M(:, :, k) = S * S' / (T-k+1);
end
A = cumsum(M, 3);
LamA = NaN(p, k0+1);
LamM = NaN(p, k0+1);

for j = 1 : k0+1
    LamA(:, j) = sort(eig(A(:, :, j)), 'descend');
    LamM(:, j) = sort(eig(M(:, :, j)), 'descend');
end
plot([0:k0]', LamA(1, :)', '*-', [0:k0]', LamM(1, :)', 'x-', [0:k0]', ...
     cumsum(LamM(1, :))', 'o-', 'Linewidth', 1);
h = legend('$\lambda_{(1)}$ of $\sum_{i=0}^k C_i C_i^T$', ...
           '$\lambda_{(1)}$ of $C_k C_k^T$',...
           '$\sum_{i=0}^k \lambda_{(1)}^{(i)}$');
set(h, 'Interpreter', 'latex', 'Fontsize', 14, 'Location', 'Southeast');
xlabel('k');
title(name);
grid on

% ratios = lambda(2:end, :) ./ lambda(1:end-1, :);
% plot(1:p-1, ratios, '+-');
% % legend('k=0,...,0', 'k=0,...,1', 'k=0,...,2', 'k=0,...,3', 'k=0,...,4', ...
% %        'k=0,...,5', 'Location', 'Southeast');
% % legend('k=0,...,5', 'k=1,...,5', 'k=2,...,5', 'Location', 'Southeast');
% xlabel('i');
% ylabel('\lambda_{(i+1)}/\lambda_{(i)}');
% title(name);
% xlim([0, 5]);



