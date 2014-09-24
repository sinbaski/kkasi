clear all
close all
q = 0.1;
ends_th = NaN(15, 2);
ends_sim = NaN(15, 2);
c = 1;
for sig = [0.1, 0.2]
    q = sig;
    for phi = [0, 0.25, 0.5, 0.7070, 0.8160]
        v = sig^2 / (1 - phi^2);
        load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
                      'sig%.4f-phi%.4f.mat'], q, sig, phi));
        ev = reshape(ev, 1, prod(size(ev)));
        ends_sim(c, :) = [min(ev), max(ev)];
        z = linspace(ends_sim(c, 1), ends_sim(c, 2), 1000);
        G = LognormalGreen(z, v, q);
        [a, b] = LognormalEnds(v, q, G);
        ends_th(c, :) = [a, b];
        c = c + 1;
    end    
end



