function index = locateit(X, pos, l, func, tdiff)
if func(X, pos, pos, tdiff) == 0
    index = pos;
    return;
else
    res = func(X, pos, size(X, 1), tdiff);
    if res == 0
        index = size(X, 1);
    elseif res < 0
        index = NaN;
    else
        index = locateit_recursion(X, pos, l, func, [pos, tdiff], ...
                                   0, [NaN, NaN]);
    end
end



