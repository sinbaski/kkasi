%% Infer innovations of a MA process
function y = ma_infer(w, theta, Theta, s)
n = length(w);
done = 0;

k = 1;
phi = NaN(1, 10);
maxj = floor(log(2.0e-2) / log(abs(Theta)));

while 1
    m = min(floor(k/s), maxj);
    j = 0:m;
    coef = theta.^(k - s.*j) .* Theta.^j;
    phi(k) = sum(coef);
    
    if m == maxj
        break;
    end
    k = k + 1;
end
phi = phi(1:k);

y = NaN(n, 1);
w = [zeros(k, 1); w];
for t = n:-1:1
    y(t) = w(t + k) + phi * w([t-1:-1:t-k] + k);
end

