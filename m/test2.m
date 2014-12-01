clear all
close all
fold = 4000;
N = 600;
T = 600;
q = N/T;

ev1 = NaN(N, fold);
for sig = [0.1, 0.2, 0.5]
    for phi = [0, 0.25, 0.5, 0.7070, 0.8160]
        for a = 1:fold
            s = [randn()*sig, NaN(1, T-1)];
            for t = 2:T
                s(t) = s(t-1) * phi + randn()*sig;
            end
            R = randn(N, T)*diag(exp(s));
            C = R*R'/T;
            ev1(:, a) = eig(C);
            % dgnl1(:, a) = C(logical(eye(N)));
        end

        if (exist(sprintf(...
            '../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-sig%.4f-phi%.4f.mat', ...
            q, sig, phi)) == 2)
            load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/' ...
                          'Eig-sig%.4f-phi%.4f.mat'], q, sig, phi), 'ev');
        else
            ev = [];
        end

        % if (exist(sprintf(...
        %     '../matfys/data/sv/normal_ret/lognormal_vol/Diag-sig%.4f-phi%.4f.mat', ...
        %     sig, phi)) == 2)
        %     load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/' ...
        %                   'Diag-sig%.4f-phi%.4f.mat'], sig, phi), 'dgnl');
        % else
        %     dgnl = [];
        % end

        ev = [ev, ev1];
        % dgnl = [dgnl, dgnl1];

        save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
                      'sig%.4f-phi%.4f.mat'], q, sig, phi));
        % save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/Diag-' ...
        %               'sig%.4f-phi%.4f.mat'], sig, phi));
    end
end

quit
