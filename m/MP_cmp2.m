clear all
close all

params =  [
    0.1068, 0.8923, 2.17
    0.1017, 0.8964, 2.40
    0.1041, 0.8908, 3.00
    0.1012, 0.8885, 4.00
          ];

spec = cellstr(['b  '; 'c  '; 'g  '; 'm  '; 'r  '; 'k  ';...
                'b--'; 'c--'; 'g--'; 'm--'; 'r--'; 'k--';...
                'b: '; 'c: '; 'g: '; 'm: '; 'r: '; 'k: ';...
                'b-.'; 'c-.'; 'g-.'; 'm-.'; 'r-.'; 'k-.';...
                ]);
q = 50;
left_end = (1 - sqrt(1/q))^2;
right_end = (1 + sqrt(1/q))^2;
x2 = linspace(left_end, right_end, 1000);
y2 = MarcenkoPasturPDF(x2, [1/q, 1]);

for n = 1 : size(params, 1)
    subplot(2, 2, n);
    texts= {};
    hold on
    c = 1;
    for N = [10, 250]
        load(sprintf(['../data/tail_exponent_%.2f/' ...
                      'GarchWishartOrdN%dQ%dA%.4fB%.4fEig.mat'],...
                     params(n, 3), N, q, params(n, 1), ...
                     params(n, 2)), 'ev');
        ev = reshape(ev, prod(size(ev)), 1);
        [x1, y1] = epdf(ev, 4, min(ev), max(ev), 300, spec{c});
        texts{c} = sprintf('#series = %d', N);
        c = c + 1;
    end
    plot(x2, y2, 'k-', 'LineWidth', 2);
    texts{c} = 'MP';
    title(sprintf('\\alpha = %.2f', params(n, 3)));
    hold off
    legend(texts);
end
