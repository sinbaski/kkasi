% interpretation of res:
% negative: index is left to the target
% 0: index is at the target
% positive: index is right to the target
function res = intra_time_sort(tid, pos, index, tdiff)
n = etime(tid(index, :), tid(pos, :));
if n < tdiff
    res = -1;
elseif n == tdiff
    res = 0;
elseif etime(tid(index-1, :), tid(pos, :)) < tdiff
    res = 0;
else
    res = 1;
end
