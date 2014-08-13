clear all
close all

N = 250;
q = 24;
T = N * q;
a = 0.94;

load(sprintf(['../data/IGarchWishartN%dQ%dA%.4fEig.mat'], ...
             N, q, a), 'ev');

ev = reshape(ev, 1, prod(size(ev)));
[x1, y1] = epdf(ev, 4, min(ev), max(ev), 200, '');
plot(x1, y1);
