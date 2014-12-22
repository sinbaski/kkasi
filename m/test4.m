clear all
close all

phi = [0
       0.25
       0.5
       0.7070
       0.8160
      ];
sig = [0.1, 0.2, 0.5];
spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
                'b: '; 'c: '; 'g: '; 'm: '; 'r: '; 'k: ';...
                'b-.'; 'c-.'; 'g-.'; 'm-.'; 'r-.'; 'k-.';...
                ]);
texts = {};
c = 1;
% hold on
lam1 = NaN(length(sig), length(phi), 2);
lam2 = NaN(length(sig), length(phi), 2);
variance = NaN(length(sig), length(phi));

q = 0.1;
for m = 1 : length(sig)
    for n = 1 : length(phi)
        % for sig = 0.1
        v = sig(m)^2/(1 - phi(n)^2);
        if v >= 0.1
            continue;
        end
        variance(m, n) = v;
        load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
                      'sig%.4f-phi%.4f.mat'], q, sig(m), phi(n)), 'ev');

        % ev = reshape(ev, 1, prod(size(ev)));
        eigmax = max(ev)';
        eigmin = min(ev)';
        lambda_max = (1+sqrt(q))^2 * (1 + 2*v*(2*sqrt(q)+1));
        lambda_min = (1-sqrt(q))^2 * (1 - 2*v*(2*sqrt(q)-1));
        subplot(3, 3, c);
        qqplot(eigmax);
        grid on
        title(sprintf(['average=%.2e v=%.4f\n'...
                      ' mean=%.2e, std=%.2e'], ...
                      mean(eigmax), v, lambda_max, std(eigmax)));
        c = c + 1;
        % lam2(m, n, :) = [mean(eigmax), (1+sqrt(q))^2 * ...
        %                  (1 + 2*v*(2*sqrt(q)+1))];
        % lam1(m, n, :) = [mean(eigmin), (1-sqrt(q))^2 * ...
        %                  (1 - 2*v*(2*sqrt(q)-1))];
        
    end
end

