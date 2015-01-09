clear all
n = 1000;
U = rand(n, 1);
X = sort(tan(pi * (U - 1/2)), 'descend');
Hill = NaN(1, n/4);
Pickands = NaN(1, n/4);
for k = 1:n/4
    Hill(k) = mean(log(X(1:k) ./ X(k)));
    Pickands(k) = log((X(k) - X(2*k))/(X(2*k) - X(4*k)))/log(2);
end
K = [1:n/4];
plot(K, Hill, K, Pickands, K, ones(1, length(K)), 'LineWidth', 2);
legend('Hill', 'Pickands', '\xi = 1');
xlabel('k: # upper order statistics');
ylabel('Estimator');
grid on
[var(Hill), var(Pickands)]
