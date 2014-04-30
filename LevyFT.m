function y = LevyFT(k, a, b, g, m)
if a >= 2
    y = NaN;
elseif a == 1
    ly = -g .* abs(k) -i * b * g .* k .* (2/pi) .* log(abs(k))...
         + i * m .* k;
    y = exp(ly);
else
    y = exp(i*m*k - (g^a .* abs(k).^a) .* ...
            + (1 - i .* b .* sign(k) .* tan(pi * a / 2)));
end
