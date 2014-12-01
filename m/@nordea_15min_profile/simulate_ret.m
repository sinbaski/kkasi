function [ret, lv] = simulate_ret(me, N, np)
%% Simulate np paths of realizations of the model.
% Each path has length N
s = me.s;
z = randn(N+(s+1), np);
y = me.jsp.Xi + me.jsp.Lambda * ...
    sinh((z - me.jsp.Gamma)./me.jsp.Delta);
w = NaN(N, np);
lv = [me.lv(end-s:end); NaN(N, np)];
for t = s+2:N+s+1
    w(t-s-1, :) = y(t, :) - me.MA * y(t - [1:s+1], :);
end
% An outer-product to augment the vector into a
% matrix of np identical columns.
for t = s+2:N+s+1
    lv(t, :) = w(t-s-1, :) + 0.99*lv(t-1, :) ...
        + 0.99*lv(t-s, :) - 0.99^2*lv(t-s-1, :);
end
lv = lv(s+2:end, :);

lv = (lv - mean(lv)) ./ std(lv) .* std(me.lv) + mean(me.lv);
ret = randn(N, np) .* exp(lv);
