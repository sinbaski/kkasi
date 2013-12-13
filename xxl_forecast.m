function w = xxl_forecast(ma, y, h)
%% Forecast h horizons into the future of x
% ma: moving average parameters
% x:  time series
% h:  number of horizons
w = NaN(h, 1);
for j = 1:h
    w(j) = ma * expectation(y, length(y)+j-[1:length(ma)]);
end
