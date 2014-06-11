clear all
close all

fold = 1e+4;
AR = [0:0.1:0.9];
tau = -log(2)./log(AR);

N = 50;
% rho = 0.5;
% q = 12;

% N = 2;
% theta = pi/4;
% q = 3000;

% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.15, 0.84], 'distr', ...
%               'Gaussian', 'TailExponent', 2.9664);
dist = struct('name', 'Garch1_1', 'prmt', [0.01, 0.15, 0.84], 'distr', ...
              'Gaussian', 'TailExponent', 2.9664);

spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
                'b: '; 'c: '; 'g: '; 'm: '; 'r: '; 'k: ';...
                'b-.'; 'c-.'; 'g-.'; 'm-.'; 'r-.'; 'k-.';...
                ]);
hold on
c = 1;
% for phi = 0
% for phi = [0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.9]
for phi = [0]
    %% The matrix elements 
    % load(sprintf('GarchWishartProperN%dQ%dRho%.3fC-%.3f.mat', ...
    %              N, q, rho, phi), 'C');
    % x = NaN(1, 2*(N-1)*fold);
    % for k = 1:fold
    %     x((k-1)*2*(N-1)+1 : k*2*(N-1)) = ...
    %         [C(1, 2:end, k), C(2:end, 1, k)'] ./ C(1, 1, k);
    % end
    % [u, v] = epdf(x, 4, 0, 1, 320, spec{c});
    % param = stblfit(x(1:1e+5), 'ecf');
    % v1 = stblpdf(u, param(1), param(2), param(3), param(4));
    % plot(u, log10(v), u, log10(v1));

    % load(sprintf('GarchWishartProperN%dQ%dRho%.3fC-%.3f.mat', ...
    %              N, q, rho, phi), 'C');
    % % C11 = reshape(C(1, 1, :), 1, fold);
    % C12 = reshape(C(1, 2, :), 1, size(C, 3));
    % C21 = reshape(C(2, 1, :), 1, size(C, 3));
    % % C22 = reshape(C(2, 2, :), 1, fold);
    
    % A = [C12, C21];
    % fprintf('var = %.3e, skew = %.2f\n', var(A), skewness(A));
    % % epdf(C11, 2, 1., -0.496, 160, spec{c});
    % % m = mean(A);
    % if c == 1
    %     subplot(1, 2, 1);
    %     grid on
    %     hold on;
    %     xlabel('C_{12} and C_{21}');
    %     ylabel('probability density function');
    % elseif c == 6
    %     legend(lgd);
    %     subplot(1, 2, 2);
    %     grid on
    %     hold on;
    %     xlabel('C_{12} and C_{21}');
    %     ylabel('probability density function');
    % end
    % if c < 6
    %     [x, y] = epdf(A, 4, -0.505, -0.495, 200, spec{c});
    % else
    %     [x, y] = epdf(A, 4, -0.505, -0.495, 200, spec{c-5});
    % end
    % [y, x] = ecdf(A);
    % I = x > -0.5;
    % stairs(log10(x(I) + 0.5), log10(1 - y(I)), spec{c});

    
    %% The eigenvalues 
    % load(sprintf('GaussianWishartN%dQ%dRho%.3fEig-%.3f.mat', ...
    %              N, q, rho, phi), 'ev');
    % E1 = ev(4, :);
    % m = mean(E1);
    % sig = std(E1);
    % fprintf('std = %.3e, skew = %.2f\n', sig, skewness(E1));
    % [x, y] = epdf(E1, 4, m-5*sig, m+5*sig, 300, spec{c});

    load(sprintf('GarchWishartProperN%dQ%dRho%.3fC-%.3f.mat', ...
                 N, q, rho, phi), 'C');
    E = [...
        reshape(C(1, 1, :), 1, size(C, 3));
        reshape(C(2, 2, :), 1, size(C, 3));
        reshape(C(3, 3, :), 1, size(C, 3));
        reshape(C(3, 3, :), 1, size(C, 3));
        ];
    % E1 = ev(1, :);
    % m = mean(E1)
    % sig = std(E1)
    % fprintf('std = %.3e, skew = %.2f\n', sig, skewness(E1));
    % [x, y] = epdf(E1, 4, m-3*sig, m+3*sig, 300, spec{c});
    % E = reshape(ev(1, :), 1, prod(size(ev, 2)));
    epdf(E(1, :), 2, 0, 60, 200, spec{c});
    epdf(E(2, :), 2, 0, 60, 200, spec{c+1});
    epdf(E(3, :), 2, 0, 60, 200, spec{c+2});
    epdf(E(4, :), 2, 0, 60, 200, spec{c+3});
    
    % hold on
    % load(sprintf('GarchWishartProperN%dQ%dRho%.3fC-%.3f.mat', ...
    %              N, q, rho, phi), 'C');
    % D = reshape(C(1, 1, :), [1, size(C, 3)]);
    % m = mean(D)
    % sig = std(D)
    % fprintf('std = %.3e, skew = %.2f\n', sig, skewness(D));
    % [x, y] = epdf(D, 4, m-3*sig, m+3*sig, 300, spec{c+1});

    % lgd{c} = 'Gaussian';
    % lgd{c+1} = 'Garch';
    % x = NaN(1, 2*(N-1)*fold);
    % for k = 1:fold
    %     x((k-1)*2*(N-1)+1 : k*2*(N-1)) = ...
    %         [C(1, 2:end, k), C(2:end, 1, k)'] ./ C(1, 1, k);
    % end

    % y1 = MarcenkoPasturPDF(x, [1/q, 1]);
    % plot(x, y1, 'r');
    % [y, x] = ecdf(E1);
    % I = x > 0.996 & x < 0.999;
    % plot(log(-x(I) + m), log(y(I)));
    % y1 = stblpdf(x, param(1), param(2), param(3), param(4));
    % z1 = stblcdf(x, param(1), param(2), param(3), param(4));
    % loglog(x, y, x, z1);
    % extract_xpnt(E1 - 1, 0.001, 0.003);
    % title(sprintf('skewness = %.3f', skewness(E1)));
    % lgd{c} = sprintf('\\phi = %.1f', phi);
    % lgd{c} = sprintf('\\phi = %.1f, skew = %.2f', phi, ...
                     % skewness(A));
    c = c + 2;
end
hold off
% legend(lgd);
