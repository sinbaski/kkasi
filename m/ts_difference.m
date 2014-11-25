%% y = (1 - B^s)^d x
function y = ts_difference(x, s, d)
if d == 0
    y = x;
    return;
elseif d < 0
    fprintf(['Error using ts_difference: degree of differencing must ' ...
             'be positive\n']);
    y = NaN;
    return;
end

y = x;
for n = 1:d;
    y = y(s+1:end) - y(1:end-s);
end
