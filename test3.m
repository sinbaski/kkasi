clear all
% filename = 'data/volvo_b.csv';
% price = dlmread(filename, ',', [1, 4, 3616, 4]);
% price = flipud(price); 

% ret = price2ret(price);
% r = (ret - mean(ret))/std(ret);

% x = 1.0:0.02:2.4;
% y = x;
% n = length(x);

% z = NaN(n, n);
% for u = 1:n
%     for v = 1:n
%         z(u, v) = NormalPowerLaw_pdf([x(u), y(v)], r);
%     end
% end
[x, y] = meshgrid([-5:0.2:5], [-5:0.2:5]);
z = exp(-1/2*x.^2.*exp(-2*y));
% z = NaN(length(x), length(y));
% z = sin(sqrt(x.^2 + y.^2));
% for i = 1:size(z, 1)
%     for j = 1:size(z, 2)
%         % z(i,j) = exp(-1/2*x(i)^2*exp(-2*y(j)));
%         z(i,j) = sin(sqrt(x(i)^2 + y(j)^2));
%     end
% end
% figure;
surf(x, y, z, 'FaceColor', 'green', 'EdgeColor', 'none');
camlight right; lighting phong;
xlabel('x')
ylabel('v');
zlabel('exp(-x^2 e^{-2v} / 2)')

% u1 = sqrt(2*pi) * exp(xi1.^2/2) .* xi1 .* normcdf(-xi1);
% u1 = xi1 ./ normcdf(-xi1) .* exp(-xi1.^2 / 2) / sqrt(2*pi) + 1;
