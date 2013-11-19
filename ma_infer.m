%% Infer innovations of a MA process
function y = ma_infer(w, theta, Theta, s)
n = length(w);
done = 0;

k = 1;
phi = NaN(1, 10);

while 1
    m = floor(k/s);
    j = 0:m;
    phi(k) = sum(theta.^(k - s.*j) .* Theta.^j);
    if phi(k) < 2.0e-2
        break;
    end
    k = k + 1;
end
phi = phi(1:k-1);

y = NaN(n-k+1, 1);
for t = n:-1:k
    y(t-k+1) = w(t) + phi * w(t-1:-1:t-k+1);
end

