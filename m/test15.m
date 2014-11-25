clear all
close all
q = 0.2;
ends_th = NaN(10);
ends_sim = NaN(10);
c = 1;
spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
                'b: '; 'c: '; 'g: '; 'm: '; 'r: '; 'k: ';...
                'b-.'; 'c-.'; 'g-.'; 'm-.'; 'r-.'; 'k-.';...
                ]);
hold on
coef = [];
for sig = [0.2]
    for phi = [0, 0.25, 0.5, 0.7070, 0.8160]
        v = sig.^2 ./ (1 - phi.^2);
        load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
                      'sig%.4f-phi%.4f.mat'], q, sig, phi));
        ev = reshape(ev, 1, prod(size(ev)));
        [y1, x1] = ecdf(ev);
        I = length(x1)-1000+1 : length(x1)-1;
        
        [u, lam2] = LognormalEnds(v, q);
        % y2 = (lam2 - x1(I)).^(3/2);
        plot(log10(lam2 - x1(I)), log10(1-y1(I)), spec{c});
        coef(c, :) = polyfit(log10(lam2 - x1(I)), log10(1-y1(I)), 1);
        
        % z = linspace(ends_sim(c, 1), ends_sim(c, 2), 1000);
        % G = LognormalGreen(z, v, q);
        % [a, b] = LognormalEnds(v, q);
        % ends_th(c, :) = [a, b];
        c = c + 1;
    end    
end
hold off


