function type = cmp_johnson_su(r);

q1 = [-2.1, -0.7, 0.7, 2.1];
probs = cdf('Normal', q1, 0, 1);
q2 = quantile(r, probs);
[R1, type] = johnsrnd([q1; q2], 1e5, 1);
[F1, X1] = ecdf(R1);

[f, x] = ecdf(r);
plot(x(x > -1.3), log(f(x > -1.3)), 'bx');
hold on
plot(X1, log(F1), 'r');
%stairs(X2, F2, 'g');
grid on
hold off;
legend('Empirical CDF', 'Fitting Johnson Su Dist.', 'Location', ...
       'SouthEast');
xlabel('x');
ylabel('CDF of y_t in semi-log scale, ln(P(y_t < x))');
