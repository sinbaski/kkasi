function grad = johnson_su_grad(param, y, z)
grad = NaN(length(y), length(param));
a = param(1);
b = param(2);
c = param(3);
m = param(4);

% z = a + b*asinh((y - m)./c);
% dz / da
grad(:, 1) = 1;

% dz / db
grad(:, 2) = (z - a)./b;

u = sqrt(1 + (y - m).^2 ./ c^2);

% dz / dm
grad(:, 4) = -b/c ./ u;

% dz / dc
% grad(3) = -b/c^2 .* (y - m) ./ u;
grad(:, 3) = grad(:, 4) ./ c .* (y - m);

