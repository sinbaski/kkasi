clear all
clc

genesis = 0;
filename1 = '../data/corr_eig_0.4.txt';
filename2 = '../data/corr_eig_0.6.txt';
filename = '../data/corr_eig_0.2.txt';

sigma = [
    1    0    0    0    0    0    0    0    0    0    0    0
    0.2  1    0    0    0    0    0    0    0    0    0    0
    0.2  0.2  1    0    0    0    0    0    0    0    0    0
    0.2  0.2  0.2  1    0    0    0    0    0    0    0    0
    0.2  0.2  0.2  0.2  1    0    0    0    0    0    0    0
    0.2  0.2  0.2  0.2  0.2  1    0    0    0    0    0    0
    0.2  0.2  0.2  0.2  0.2  0.2  1    0    0    0    0    0
    0.2  0.2  0.2  0.2  0.2  0.2  0.2  1    0    0    0    0
    0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  1    0    0    0
    0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  1    0    0
    0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  1    0
    0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  1
        ];
sigma = sigma + sigma' - eye(size(sigma, 1));
N = size(sigma, 1);

if genesis
    filename = '../data/corr_eig_0.2.txt';
    %N = 800;
    T = 800;
    %Q = T/N;

    for k = 1:400
        random_matrix_eig(T, sigma, filename);
    end
else
    lambda = [load(filename), load(filename1), load(filename2)];

    a = min(lambda);
    b = max(lambda);

    cen = linspace(min(a), max(b), 300);
    dx = cen(2) - cen(1);
    pdf = double(hist(lambda, cen)) ./ (N*400*dx);
    plot(cen, pdf);
end

% a1 = 1+1/Q-2*sqrt(1/Q);
% b1 = 1+1/Q+2*sqrt(1/Q);

% x = linspace(a1, b1, 200);
% f = @(u) Q./(2*pi.*u) .* sqrt((b1-u).*(u-a1));
% c = integral(f, a1, b1);

% y = Q./(2*pi.*x) .* sqrt((b1-x).*(x-a1));
% plot(x, y, 'r');
% hold on

%plot(cen, pdf(:, 1), cen, pdf(:, 2), cen, pdf(:, 3));
