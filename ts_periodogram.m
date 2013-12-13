function I = ts_periodogram(f, a)
n = length(a);
t = [1:n]';

I = NaN(length(f), 1);
for k = 1:length(f)
    I(k) = (cos(2*pi*f(k)*t)'*a)^2 + (sin(2*pi*f(k)*t)'*a)^2;
    I(k) = I(k)*2/n;
end

