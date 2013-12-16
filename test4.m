clear all
close all
f1 = @(v, xi) exp(-(v+1).^2/2 - xi^2 .* exp(-2*v) ./ 2 + 1/2) ./ ...
     (2*pi);
f2 = @(xi) integral(@(v) f1(v, xi), -Inf, Inf);

xi = [-3:0.1:3]';
y = NaN(length(xi), 1);
for i = 1:length(xi)
    if xi(i) == 0
        y(i) = exp(1/2) / sqrt(2*pi);
    else
        y(i) = f2(xi(i));
    end
end
plot(xi, y, xi, normpdf(xi, 0, 1));
