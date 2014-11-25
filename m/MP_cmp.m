clear all
close all
% data = struct('alpha', {2.17; 2.40; 3.00; 4.00; 5.00}, ...
%               'params', {
%                   [
%                       0.1068,  0.8923, NaN
%                       0.2009, 0.7962, NaN
%                       0.2991,  0.6946, NaN
%                       0.4001, 0.5889, NaN
%                       0.4996,  0.4835, NaN
%                       0.5995, 0.3786, NaN
%                       0.6992,  0.2717, NaN
%                       0.8009, 0.1608, NaN
%                       0.8970,  0.0541, NaN
%                   ]
                  
%                   [
%                       0.0663, 0.9329, NaN
%                       0.1017, 0.8964, NaN
%                       0.2016, 0.7914, NaN
%                       0.3013, 0.6838, NaN
%                       0.4007, 0.5735, NaN
%                       0.5011, 0.4593, NaN
%                       0.5992, 0.3443, NaN
%                       0.6992, 0.2221, NaN
%                       0.8013, 0.0932, NaN
%                       0.8684, 0.0040, NaN
%                   ]
                  
%                   [                  
%                       0.1041, 0.8908, NaN
%                       0.1994, 0.7822, NaN
%                       0.3008, 0.6578, NaN
%                       0.4000, 0.5273, NaN
%                       0.5000, 0.3874, NaN
%                       0.6031, 0.2259, NaN
%                       0.7156, 0.0289, NaN
%                   ]
                  
%                   [
%                       0.1012, 0.8885, NaN
%                       0.1998, 0.7592, NaN
%                       0.2996, 0.6060, NaN
%                       0.4001, 0.4239, NaN
%                       0.5004, 0.2082, NaN
%                       0.5491, 0.0784, NaN
%                   ]

%                   [
%                       0.1042, 0.8784, NaN
%                       0.2008, 0.7302, NaN
%                       0.2999, 0.5332, NaN
%                       0.3999, 0.2681, NaN
%                       0.4497, 0.1009, NaN
%                   ]
%                    });
data = struct('a1', {0.10; 0.20; 0.30},...
              'params', {
                  [
                      0.1068, 0.8923, 2.17, NaN
                      0.1017, 0.8964, 2.40, NaN
                      0.1041, 0.8908, 3.00, NaN
                      0.1012, 0.8885, 4.00, NaN
                      0.1042, 0.8784, 5.00, NaN
                  ]
                  
                  [
                      0.2009, 0.7962, 2.17, NaN
                      0.2016, 0.7914, 2.40, NaN
                      0.1994, 0.7822, 3.00, NaN
                      0.1998, 0.7592, 4.00, NaN
                      0.2008, 0.7302, 5.00, NaN
                  ]
                  
                  [
                      0.2991, 0.6946, 2.17, NaN
                      0.3013, 0.6838, 2.40, NaN
                      0.3008, 0.6578, 3.00, NaN
                      0.2996, 0.6060, 4.00, NaN
                      0.2999, 0.5332, 5.00, NaN
                  ]
                   });

spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
                'b: '; 'c: '; 'g: '; 'm: '; 'r: '; 'k: ';...
                'b-.'; 'c-.'; 'g-.'; 'm-.'; 'r-.'; 'k-.';...
                ]);

N = 50;
q = 24;
for n = 3
    figure;
    texts= {};
    hold on
    % left_end = (1 - sqrt(1/q))^2;
    % right_end = (1 + sqrt(1/q))^2;
    % x = linspace(left_end, right_end, 1000);
    % y = MarcenkoPasturPDF(x, [1/q, 1]);
    % plot(x, y, 'k-', 'LineWidth', 2);
    % texts{1} = 'MP';
    for u = 1:size(data(n).params, 1)
        % load(sprintf('../data/tail_exponent_%.2f/GarchWishartOrdN%dQ%dA%.4fB%.4fEig.mat',...
        %              data(n).alpha, N, q, data(n).params(u, 1), data(n).params(u, 2)), 'ev');
        load(sprintf('../data/tail_exponent_%.2f/GarchWishartOrdN%dQ%dA%.4fB%.4fEig.mat',...
                     data(n).params(u, 3), N, q, data(n).params(u, ...
                                                          1), data(n).params(u, 2)), 'ev');
        eigmin = min(ev);
        [x1, y1] = epdf(eigmin, 4, min(eigmin), max(eigmin), 200, '');
        % ev = reshape(ev, prod(size(ev)), 1);
        % [x1, y1] = epdf(ev, 4, min(ev), max(ev), 400, '');
        % y2 = MarcenkoPasturPDF(x1, [1/q, 1]);
        stairs(x1, y1, spec{u});
        texts{u} = sprintf('\\alpha = %.2f', data(n).params(u, 3));
        % data(n).params(u, end) = sqrt(mean((y1 - y2).^2));
        % data(n).params(u, 3) = HellingerDistance(...
        %     @(x) y1, ...
        %     @(x) MarcenkoPasturPDF(x, [1/q, 1]), x1);
    end
    hold off
    legend(texts);
end
xlim([0.52, 0.72]);
ylim([0, 40]);
grid on
xlabel('Smallest eigenvalue', 'Fontsize', 14);
ylabel('Probability density', 'Fontsize', 14);
title(sprintf('\\alpha_1 = %.1f', data(n).a1), 'Fontsize', 14);

 

% hold on
% colors = {'xb-', 'xg-', 'xc-', 'xm-', 'xr-'};
% for n = 1 : size(data, 1)
%     m = size(data(n).params, 1);
%     plot(data(n).params(:, 1), data(n).params(:, 3), colors{n});
% end
% hold off
% xlabel('\alpha_1', 'fontsize', 14);
% ylabel('Distance to MP (mean square error)', 'fontsize', 14);
% legend('\alpha = 2.17', '\alpha = 2.40', '\alpha = 3.00', '\alpha = 4.00', '\alpha = 5.00');
