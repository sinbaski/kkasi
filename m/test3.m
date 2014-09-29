clear all
close all

N = 100;
T = 1000;
q = N/T;
sig = 1;
phi = 0.707;
spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
                'b: '; 'c: '; 'g: '; 'm: '; 'r: '; 'k: ';...
                'b-.'; 'c-.'; 'g-.'; 'm-.'; 'r-.'; 'k-.';...
                ]);
texts = {};
c = 1;
hold on
q = 0.1;
%for sig = [0.1, 0.5, 1]
for sig = 0.5
    load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/q%.1f/Eig-' ...
                  'sig%.4f-phi%.4f.mat'], q, sig, 0));

    % ev = reshape(ev, 1, prod(size(ev)));
    eigmax = max(ev)';
    [lam, den] = epdf(eigmax, 1, min(eigmax), max(eigmax), 400, spec{c});
    
    % [Z2, G2] = epdf(ev, 1, min(ev), 10, 100, '');

    % lammax1= max(ev);
    [u, lammax2] = LognormalEnds(sig^2, q);
    a = LognormalEndsConstA(sig^2, q);
    lammax3 = LognormalEndsFun3(sig^2, q, a);
    % a = LognormalEndsConstA(sig^2, q);
    texts{c} = sprintf('\\lambda_{max} = %.2f', lammax2);
    % Z1 = linspace(max(min(ev),0.05), max(ev), 1000);

    
    % Z1 = linspace(lammax1, lammax2, 1000);
    % G1 = LognormalGreen(Z1, sig.^2, q);
    % plot(Z1, real(G1), spec{c});
    % texts{c} = sprintf('a = %e', a);

    % plot(Z1, -imag(G1)./pi, Z2, G2, 'LineWidth', 2);
    % stairs(Z2, G2, spec{c}, 'LineWidth', 1);
    % plot(Z1, -imag(G1)./pi, spec{c}, 'LineWidth', 1);

    % legend('Theoretical', 'Simulation');
    % texts{2*c-1} = sprintf('v = %.2f, em', sig^2);
    % texts{2*c} = sprintf('v = %.2f', sig^2);
    c = c + 1;
end

% a = (1 - sqrt(q))^2;
% b = (1 + sqrt(q))^2;
% Z1 = linspace(a, b, 1000);
% MPGreen = @(z) (z + q - 1 - sqrt((z-a).*(z-b)))./(2*q.*z);
% g = MPGreen(Z1);
% plot(Z1, -imag(g)./pi, 'r', 'LineWidth', 2);

% F1 = cumsum((-imag(g)./pi))*(Z1(2) - Z1(1));
% F1 = F1 ./ F1(end);
% plot(log10(Z1), log10(1-F1), 'r', 'LineWidth', 2);
hold off

% xlabel('Eigenvalues');
% ylabel('Probability Densities of Eigenvalues');
% xlim([0, 10]);
% ylim([0, 1.4]);
grid on
    
% [Z2, P2] = epdf(ev, 1, min(ev), max(ev), 500, '');
% plot((Z2), (P2), 'ro');
% hold on

% plot((Z1), (-imag(G1)./pi), 'b', 'LineWidth', 2);
% hold off

% grid on
legend(texts);
% ylabel('Tail function of eigenvalue distribution');


