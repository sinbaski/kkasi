function e = expectation(y, t)
%% The expectation of y
% t: vector or scalar. The points at which the expectation is
% evaluated.
e = zeros(length(t), 1);
n = length(y);
idx = t <= n & t >= 1;
t = t(idx);
e(idx) = y(t);
