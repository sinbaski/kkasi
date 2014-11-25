clear all
close all
dof = 3;
d = (4-dof)/(2*(dof - 1)) * 0.9;
p = 10;
n = floor(p^(1/d));
fold = 5000;


load(sprintf(['../matfys/data/student3/Eig.mat']), 'ev');
% ev = [];
ev_new = NaN(p, fold);
for c = 1 : fold
    R = trnd(dof, p, n);
    C = R * R';
    if dof > 2
        C = C - eye(p) .* (dof / (dof - 2));
    end
    ev_new(:, c) = eig(C);
end
ev = [ev, ev_new];
save(sprintf(['../matfys/data/student3/Eig.mat']), 'ev');
quit

