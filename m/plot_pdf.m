function plot_pdf(data, N, theopdf, lag, idx, loglog, spec1, spec2)
x = linspace(min(data)*0.8, max(data)*1.2, N+lag-1);
y = hist(data, x)/length(data)/(x(2)-x(1));
x = tsmovavg(x, 's', lag);
y = tsmovavg(y, 's', lag);
x = x(lag:end);
y = y(lag:end);
if ~loglog
    if ~isempty(spec1) && ~isempty(spec2)
        plot(x(idx), y(idx), spec1, x(idx), theopdf(x(idx)), spec2);
    else
        plot(x(idx), y(idx), x(idx), theopdf(x(idx)));
    end
else
    if ~isempty(spec1) && ~isempty(spec2)
        plot(log(x(idx)), log(y(idx)), spec1, log(x(idx)), log(theopdf(x(idx))), spec2);
    else
        plot(log(x(idx)), log(y(idx)), log(x(idx)), log(theopdf(x(idx))));
    end
end
