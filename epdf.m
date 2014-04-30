function [x, y] = epdf(data, lag, x1, x2, N, spec)
% x1 = min(min(data1), min(data2));
% x2 = max(max(data1), max(data2));
x = linspace(x1, x2, N);
y = hist(data, x) ./ length(data) ./ (x(2) - x(1));
x = x(2:end-1);
y = y(2:end-1);

x = tsmovavg(x, 's', lag);
x = x(lag:end);

y = tsmovavg(y, 's', lag);
y = y(lag:end);

if ~isempty(spec)
    plot(log10(x), log10(y), spec, 'LineWidth', 2);
end
