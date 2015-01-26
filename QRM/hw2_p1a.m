clear all
close all
mysql = get_mysql();
X = NaN(1700, 2);
assets ={'DAX', 'BMW_DE'};
for k = 1 : size(X, 2)
    data = cell2mat(fetch(mysql, sprintf(['select closing from %s order ' ...
                        'by day desc limit 1701'], assets{k})));
    X(:, k) = price2ret(flipud(data));
end
close(mysql);
C = cov(X);
B = chol(C)';
Y = inv(B) * (X - repmat(mean(X), size(X, 1), 1))';

T1 = trnd(3.5, 2, size(Y, 2));
X1 = B * T1 + repmat(mean(X)', 1, size(T1, 2));
L1 = -[100, 100] * X1;


% % angles = NaN(1, size(Z, 1));
% R = NaN(1, size(Z, 1));
% for k = 1 : size(Z, 2)
%     R(k) = norm(Z(:, k));
%     % Z(:, k) = Z(:, k) ./ norm(Z(:, k));
%     % angles(k) = angle(Z(1, k) + i*Z(2, k));
% end
