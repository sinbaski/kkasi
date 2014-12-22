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

for q = [0.1, 0.2, 0.5, 1]
    for m = 1 : length(sig)
        for n = 1 : length(phi)
            % for sig = 0.1
            v = sig(m)^2/(1 - phi(n)^2);
            variance(m, n) = v;
            load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
                          'sig%.4f-phi%.4f.mat'], q, sig(m), phi(n)), 'ev');

            % ev = reshape(ev, 1, prod(size(ev)));
            eigmax = max(ev)';
            eigmin = min(ev)';
            
            lam2(m, n, :) = [max(eigmax), (1+sqrt(q))^2 * ...
                             (1 + 2*v*(2*sqrt(q)+1))];
            lam1(m, n, :) = [min(eigmin), (1-sqrt(q))^2 * ...
                             (1 - 2*v*(2*sqrt(q)-1))];
            
        end
    end

    [V, I] = sort(reshape(variance, numel(variance), 1));
    A = NaN(numel(lam1(:, :, 1)), 2);
    B = NaN(numel(lam2(:, :, 1)), 2);

    A(:, 1) = reshape(lam1(:, :, 1), numel(lam1(:, :, 1)), 1);
    A(:, 1) = A(I, 1);
    A(:, 2) = reshape(lam1(:, :, 2), numel(lam1(:, :, 2)), 1);
    A(:, 2) = A(I, 2);
    B(:, 1) = reshape(lam2(:, :, 1), numel(lam2(:, :, 1)), 1);
    B(:, 1) = B(I, 1);
    B(:, 2) = reshape(lam2(:, :, 2), numel(lam2(:, :, 2)), 1);
    B(:, 2) = B(I, 2);

    % f = @(param, X) param(1) .* exp(param(2) .* X) + param(3);
    % param1 = lsqcurvefit(f, [1, 1, (1 - sqrt(q))^2], V, A(:, 1));
    % param2 = lsqcurvefit(f, [1, 1, (1 + sqrt(q))^2], V, B(:, 1));

    % X = linspace(min(V), max(V), 500)';
    % semilogx(V, A(:, 1), '+', X, f(param1, X), V, A(:, ...
    %                                                   2), '-o');
    subplot(2, 2, c);
    semilogx(V, A(:, 1), 'g-', V, A(:, 2), 'r-');
    grid on
    % legend('simulated E(\lambda_{min})', 'Fit to exponential', ['Linear ' ...
    %                    'approximation'], 'Location', 'Northwest');
    xlabel('log(v)');
    ylabel('min(\lambda_{max})');
    title(sprintf('q = %.1f', q));
    c = c + 1;
end


