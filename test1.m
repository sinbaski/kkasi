M = central(q);
M(1) = f1(q);
M(2) = sqrt(M(2));
M(3) = M(3)/M(2)^3;
M(4) = M(4)/M(2)^4;
fprintf('MLE parameters:\n');
fprintf('%.4e & %.4e & %.4e & %.4e\n', q);
fprintf('Sample Moments:\n');
fprintf('%.4e & %.4e & %.4e & %.4e\n', mean(ret), std(ret), ...
        skewness(ret), kurtosis(ret));
fprintf('Theoretical Moments:\n');
fprintf('%.4e & %.4e & %.4e & %.4e\n', M);

[y, x] = ecdf(ret);
x1 = linspace(min(ret), max(ret), 500);
y1 = slv_cdf((x1 - q(4)).* exp(-q(3)), -0.0157, q(2));
% y2 = slv_cdf((x1 - p(4)).* exp(-p(3)), p(1), p(2));

subplot(1, 2, 1);
plot(log(-x(x<0)), log(y(x<0)), ...
     log(-x1(x1<0)), log(y1(x1<0)));
grid on
xlabel('ln(-x) where x < 0', 'Interpreter', 'tex');
ylabel('ln(P(r < x)) where x < 0', 'Interpreter', 'tex');
legend('empirical CDF', 'Model CDF', 'Location', 'Southwest');

subplot(1, 2, 2);
plot(log(x(x>0)), log(1-y(x>0)), log(x1(x1>0)), log(1-y1(x1>0)));
grid on
xlabel('ln(x) where x > 0', 'Interpreter', 'tex');
ylabel('ln(P(r > x)) where x > 0', 'Interpreter', 'tex');
legend('empirical compl. CDF', 'Model compl. CDF', 'Location', ...
       'Southwest');
