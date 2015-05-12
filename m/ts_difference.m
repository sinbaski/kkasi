%% y = (1 - B^s)^d x
function y = ts_difference(x, SD)
y = x;
for a = 1 : size(SD, 1)
    for b = 1 : SD(a, 2)
        y = y(1+SD(a,1):end) - y(1:end-SD(a,1));
    end
end
