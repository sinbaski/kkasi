clear all
close all

sig = [0.1, 0.2, 0.5];
spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
                'b: '; 'c: '; 'g: '; 'm: '; 'r: '; 'k: ';...
                'b-.'; 'c-.'; 'g-.'; 'm-.'; 'r-.'; 'k-.';...
                ]);
texts = {};
c = 1;
hold on
q = 0.1;
for m = 1 : length(sig)
    % for sig = 0.1
    v = sig(m)^2;
    load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
                  'sig%.4f-phi%.4f.mat'], q, sig(m), 0), 'ev');

    % ev = reshape(ev, 1, prod(size(ev)));
    % eigmax = max(ev)';
    eigmin = min(ev)';
    epdf(eigmin, 4, min(eigmin), max(eigmin), 160, spec{c});
    texts{c} = sprintf('v = %.2f', v);
    grid on;
    c = c + 1;
end
legend(texts);
