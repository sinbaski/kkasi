clear all
close all
% mysql = get_mysql();
X = NaN(1700, 2);
% assets ={'DAX', 'BMW_DE'};
% for k = 1 : size(X, 2)
%     data = cell2mat(fetch(mysql, sprintf(['select closing from %s order ' ...
%                         'by day desc limit 1701'], assets{k})));
%     X(:, k) = price2ret(flipud(data));
% end
% close(mysql);
A = importdata('DAX.txt');
X(:, 1) = A(:, 2);

A = importdata('bmw_returns.txt');
X(:, 2) = A(:, 2);

C = cov(X);
B = chol(C)';
Y = inv(B) * (X - repmat(mean(X), size(X, 1), 1))';
Y = Y';

% Z = NaN(size(Y));
% angles = NaN(size(Y, 1), 1);
% for k = 1 : size(Z, 1)
%     Z(k, :) = Y(k, :) ./ norm(Y(k, :));
%     angles(k) = angle(Z(k, 1) + i*Z(k, 2));
% end
% Tsim = trnd(3.2, size(Y, 1), 2);

Tsim = trnd(3.2, 40000, 2);
Xsim = Tsim * B';
W = 100 * ones(1, 2);
L = -Xsim * W';


% X1 = B * T1 + repmat(mean(X)', 1, size(T1, 2));
% L1 = -[100, 100] * X1;


% % angles = NaN(1, size(Z, 1));
% R = NaN(1, size(Z, 1));
% for k = 1 : size(Z, 2)
%     R(k) = norm(Z(:, k));
%     % Z(:, k) = Z(:, k) ./ norm(Z(:, k));
%     % angles(k) = angle(Z(1, k) + i*Z(2, k));
% end
