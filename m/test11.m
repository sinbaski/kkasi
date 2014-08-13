clear all
close all

fold = 3999;

N = 250;
q = 24;
T = N * q;
a = 0.94;
dist = struct('name', 'I-Garch1', 'prmt', a);


if (exist(sprintf(['../data/IGarchWishartN%dQ%dA%.4fEig.mat'], ...
                  N, q, a), 'file') == 2)
    load(sprintf(['../data/IGarchWishartN%dQ%dA%.4fEig.mat'], ...
                  N, q, a), 'ev');
else
    ev = [];
end

ev1 = NaN(N, fold);
for i = 1:fold
    R = gen_ret_mtx(N, T, dist, 0, []);
    C1(:, :, i) = R*R' ./ T;
    ev1(:, i) = eig(C1(:, :, i));
end
ev = [ev, ev1];
save(sprintf(['../data/IGarchWishartN%dQ%dA%.4fEig.mat'], ...
                  N, q, a), 'ev');
fprintf('Saved files N%d Q%d %.2f.\n', N, q, a);

