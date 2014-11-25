function C = ts_cumm_periodogram(f, a)
n = length(a);
N = 200;
dg = 1/N;
for k = 1:length(f)
    x = 0;
    for g = dg:dg:f(k)
        x = x + ts_periodogram(g, a)/n;
    end
    C(k) = x/var(a);
end
