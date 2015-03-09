clear all
close all
fold = 1.0e+4;
N = 600;
T = 600;
q = N/T;
phi = 0;

ev = NaN(fold);
eigvec = NaN(N, fold);
diagelem = NaN(fold);

% for sig = sqrt([0.02, 0.03, 0.5, 1.0])
for sig = sqrt(0.5)
    for a = 1:fold
        s = randn(1, T)*sig;
        % for t = 2:T
        %     s(t) = s(t-1) * phi + randn()*sig;
        % end
        R = randn(N, T)*diag(exp(s));
        C = R*R'/T;
        [V, D] = eig(C);
        [v, id] = max(D(logical(eye(N))));
        ev(a) = v;
        eigvec(:, a) = V(:, id);
        diagelem(a) = max(C(logical(eye(N))));
    end

    % if (exist(sprintf(...
    %     '../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-sig%.4f-phi%.4f.mat', ...
    %     q, sig, phi)) == 2)
    %     load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/' ...
    %                   'Eig-sig%.4f-phi%.4f.mat'], q, sig, phi), 'ev');
    % else
    %     ev = [];
    % end

    % if (exist(sprintf(...
    %     '../matfys/data/sv/normal_ret/lognormal_vol/Diag-sig%.4f-phi%.4f.mat', ...
    %     sig, phi)) == 2)
    %     load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/' ...
    %                   'Diag-sig%.4f-phi%.4f.mat'], sig, phi), 'dgnl');
    % else
    %     dgnl = [];
    % end

    % ev = [ev, ev1];
    % dgnl = [dgnl, dgnl1];

    save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
                  'sig%.4f-phi%.4f.mat'], q, sig, phi), 'ev');
    save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eigvec-' ...
                  'sig%.4f-phi%.4f.mat'], q, sig, phi), 'eigvec');
    save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Diagelem-' ...
                  'sig%.4f-phi%.4f.mat'], q, sig, phi), 'diagelem');
    % save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/Diag-' ...
    %               'sig%.4f-phi%.4f.mat'], sig, phi));
end

quit
