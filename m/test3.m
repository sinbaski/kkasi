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
for sig = [0.1, 0.5, 1]
    % for sig = 0.5
    % load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/Eig-' ...
    %               'sig%.4f-phi%.4f.mat'],sig, 0));

    % ev = reshape(ev, 1, prod(size(ev)));

    % F1 = cumsum((-imag(G1)./pi))*(Z1(2) - Z1(1));
    % F1 = F1 ./ F1(end);

    % [F2, Z2] = ecdf(ev);
    % plot(Z2, F2);
    
    % plot(log10(Z2(2:end-1)), log10(1- F2(2:end-1)), ...
    %      spec{c},'LineWidth', 2);
    
    % [Z2, G2] = epdf(ev, 1, min(ev), max(ev), 1000, '');
    Z1 = linspace(0.1, 10, 1000)';
    G1 = LognormalGreen(Z1, sig.^2, q);

    % plot(Z1, -imag(G1)./pi, Z2, G2, 'LineWidth', 2);
    plot(Z1, -imag(G1)./pi, spec{c}, 'LineWidth', 2);
    grid on
    % legend('Theoretical', 'Simulation');
    % xlabel('Eigenvalues');
    % ylabel('Probability Densities of Eigenvalues');
    % texts{c} = sprintf('\\phi = %.4f', phi);
    c = c + 1;
end

a = (1 - sqrt(q))^2;
b = (1 + sqrt(q))^2;
Z1 = linspace(a, b, 1000);
MPGreen = @(z) (z + q - 1 - sqrt((z-a).*(z-b)))./(2*q.*z);
g = MPGreen(Z1);
% F1 = -imag(g)./pi;
plot(Z1, -imag(g)./pi, 'r', 'LineWidth', 2);

% F1 = cumsum((-imag(g)./pi))*(Z1(2) - Z1(1));
% F1 = F1 ./ F1(end);
% plot(log10(Z1), log10(1-F1), 'r', 'LineWidth', 2);
hold off
% [Z2, P2] = epdf(ev, 1, min(ev), max(ev), 500, '');
% plot((Z2), (P2), 'ro');
% hold on

% plot((Z1), (-imag(G1)./pi), 'b', 'LineWidth', 2);
% hold off

% grid on
% legend(texts);
% ylabel('Tail function of eigenvalue distribution');


