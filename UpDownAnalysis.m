clear all
filename = 'data/OMXS30.csv';
price = dlmread(filename, ';', [1, 1, 6857, 3]);
price = flipud(price); 
high = price(:, 1);
low = price(:, 2);
closing = price(:, 3);

closing = closing(high > low & closing > 0);
ret = price2ret(closing);

s = sign(ret(1));
c = 1;
ups = zeros(1, 20);
downs = zeros(1, 20);
for i = 2:length(ret)
    t = sign(ret(i));
    if t == s || t == 0
        c = c + 1;
    else
        if s > 0
            ups(c) = ups(c) + 1;
        else
            downs(c) = downs(c) + 1;
        end
        s = t;
        c = 1;
    end
end

Pu = ups(1:max(find(ups > 0)))./sum(ups);
Pd = downs(1:max(find(downs > 0)))./sum(downs);
plot(1:length(Pu), Pu, 'g', 1:length(Pd), Pd, 'r');
grid on
