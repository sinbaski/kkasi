clear all
close all
N = 100;
T = 1000;
fold = 40000;
q = N/T;
sig = 0.5;
ev1 = NaN(N, fold);
dgnl1 = NaN(N, fold);

for phi = [0, 0.707 0.816 0.866]
    for a = 1:fold
        s = [randn(), NaN(1, T-1)]*sig;
        for t = 2:T
            s(t) = s(t-1) * phi + randn()*sig;
        end
        R = randn(N, T)*diag(exp(s));
        C = R*R'/T;
        ev1(:, a) = eig(C);
        dgnl1(:, a) = C(logical(eye(N)));
    end

    if (exist(sprintf(...
        '../matfys/data/sv/normal_ret/lognormal_vol/Eig-sig%.4f-phi%.4f.mat', ...
        sig, phi)) == 2)
        load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/' ...
                      'Eig-sig%.4f-phi%.4f.mat'], sig, phi), 'ev');
    else
        ev = [];
    end

    if (exist(sprintf(...
        '../matfys/data/sv/normal_ret/lognormal_vol/Diag-sig%.4f-phi%.4f.mat', ...
        sig, phi)) == 2)
        load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/' ...
                      'Diag-sig%.4f-phi%.4f.mat'], sig, phi), 'dgnl');
    else
        dgnl = [];
    end

    ev = [ev, ev1];
    dgnl = [dgnl, dgnl1];

    save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/Eig-' ...
                  'sig%.4f-phi%.4f.mat'], sig, phi));
    save(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/Diag-' ...
                  'sig%.4f-phi%.4f.mat'], sig, phi));
end
quit
