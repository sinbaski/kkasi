clear all
m = 902;
l = 49;

k = 1:m/2;
s1 = prod(2*k ./ (2*k-1));

k=1:(l+1)/2;
s2 = prod((m+2*k-2) ./ (m+l+2*k));

s = s1 * s2 / m / pi;
