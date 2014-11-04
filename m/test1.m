clear all
close all

mysql = get_mysql();
stmt = sprintf('select symbol from DAX_components;');
data = fetch(mysql, stmt);
p = 29;
T = 3000;

T1 = 300;
periods = T/T1;

R = NaN(p, T);
counts = NaN(1, p);
for k = 1 : p
    closing = fetch(mysql, sprintf(['select closing from %s order by day desc ' ...
                        'limit %d;'], strrep(data{k}, '.', '_'), ...
          T));
    R(k, :) = cell2mat(flipud(closing))';
end
close(mysql);

lambda = NaN(p, periods);
for k = 1: periods
    C = R(:, (k-1)*T1 + 1 : k*T1) * R(:, (k-1)*T1 + 1 : k*T1)';
    lambda(:, k) = sort(eig(C), 'descend');
end
I = lambda(1, :) > 0 & lambda(2, :) > 0;
A = log10(lambda(1, I)./lambda(2, I));
A = A./std(A);

B = log10(1./rand(1, 1000));
B = B./std(B);
qqplot(A, B);

% results = NaN(1, p-1);
% U = sort(rand(500, p-1), 'descend');
% for k = 1: p-1
%     eig_sample = log(lambda(p, :)./lambda(p-k, :));
%     eig_sample = eig_sample./std(eig_sample);
%     uni_sample = log(U(k, :));
%     uni_sample = uni_sample ./ std(uni_sample);
%     results(k) = kstest2(eig_sample, uni_sample);
% end

% fd = fopen('results.txt', 'wt');
% for k = 1 : p - 1
%     fprintf(fd, '%d, ', results(k));
%     fprintf(fd, '\n');
% end
% fclose(fd);
% Check whether lambda()
% C = R*R';
% lambda = flipud(sort(eig(C)));
% lambda = lambda ./ lambda(1);


