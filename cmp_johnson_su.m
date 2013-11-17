function type = cmp_johnson_su(r);

q1 = [-2.1, -0.7, 0.7, 2.1];
probs = cdf('Normal', q1, 0, 1);
q2 = quantile(r, probs);
[R1, type] = johnsrnd([q1; q2], 1e5, 1);
[F1, X1] = ecdf(R1);

ecdf(r);
hold on
stairs(X1, F1, 'r');
%stairs(X2, F2, 'g');
grid on
hold off;
