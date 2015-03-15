clear all
close all
fold = 1.0e+2;
N = 600;
T = 600;
q = N/T;
phi = 0;

ev = NaN(2, fold);
eigvec = NaN(N, fold, 2);
diagelem = NaN(2, fold);

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
        % [eigval, I] = sort(D(logical(eye(N))), 'descend');
        E = sort(C(logical(eye(N))));
        
        ev(1, a) = D(N, N);
        ev(2, a) = D(N/2, N/2);

        eigvec(:, a, 1) = V(:, N);
        eigvec(:, a, 2) = V(:, N/2);

        diagelem(1, a) = E(N);
        diagelem(2, a) = E(N/2);
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

    % save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
    %               'sig%.4f-phi%.4f.mat'], q, sig, phi), 'ev');
    % save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eigvec-' ...
    %               'sig%.4f-phi%.4f.mat'], q, sig, phi), 'eigvec');
    % save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Diagelem-' ...
    %               'sig%.4f-phi%.4f.mat'], q, sig, phi), 'diagelem');
end

% quit
