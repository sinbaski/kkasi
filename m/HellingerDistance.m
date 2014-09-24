% Calculate the Hellinger Distance between 2 density functions
function d = HellingerDistance(f, g, x)
x = reshape(x, length(x), 1);
dx = x(2:end) - x(1:end-1);
dx = [dx; dx(end)];
d = sqrt(reshape(f(x), size(dx')) .* ...
         reshape(g(x), size(dx'))) * dx;
d = sqrt(1 - d)*sqrt(2);



