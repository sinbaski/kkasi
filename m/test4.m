clear all
close all

a = linspace(10^(-2.7), 10^(-2.65), 8000);
for q = 0.1
    for v = 0.5
        % g = LognormalGreen(mean(eigmax), v, q);
        I = LognormalEndsFun1(v, q, a);
        plot(log10(a), log10(I), log10(a), -2*log10(a));
    end
end
xlabel('log_{10}a');
legend('the integral', 'y=1/a^2');
grid on
