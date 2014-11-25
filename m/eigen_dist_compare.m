clear all
close all

hold on;

load(sprintf('GaussianWishartEig005-%d.mat', 1), 'ev');
ev = reshape(ev, prod(size(ev)), 1);
epdf(ev, 6, min(ev)*0.9, max(ev)*1.1, 240, 'r');

load(sprintf('GarchWishartEig-%d.mat', 1), 'ev');
ev = reshape(ev, prod(size(ev)), 1);
epdf(ev, 6, 0, max(ev)*1.1, 240, 'b');


load(sprintf('FatWishartEig-%d.mat', 1), 'ev');
ev = reshape(ev, prod(size(ev)), 1);
epdf(ev, 6, min(ev)*0.9, max(ev)*1.1, 240, 'g');

hold off
