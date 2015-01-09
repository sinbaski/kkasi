clear all
X = importdata('sp500_1928-2000.txt');
X = reshape(X, numel(X), 1);
gains = sort(X(X > 0), 'descend');
losses = sort(-X(X < 0), 'descend');
A = sort(abs(X), 'descend');

Hg = NaN(floor(length(gains)/10), 1);
Hl = NaN(floor(length(losses)/10), 1);
Ha = NaN(floor(length(A)/10), 1);

Cg = NaN(floor(length(gains)/10), 1);
Cl = NaN(floor(length(losses)/10), 1);
Ca = NaN(floor(length(A)/10), 1);

Qg = NaN(floor(length(gains)/10), 1);
Ql = NaN(floor(length(losses)/10), 1);
Qa = NaN(floor(length(A)/10), 1);

for k = 1 : length(Hg)
    Hg(k) = 1/mean(log(gains(1:k)./gains(k)));
    Cg(k) = k/length(gains) * gains(k)^Hg(k);
    Qg(k) = (1.0e-2/Cg(k))^(-1/Hg(k));
end

for k = 1 : length(Hl)
    Hl(k) = 1/mean(log(losses(1:k)./losses(k)));
    Cl(k) = k/length(losses) * losses(k)^Hl(k);
    Ql(k) = (1.0e-2/Cl(k))^(-1/Hl(k));
end

for k = 1 : length(Ha)
    Ha(k) = 1/mean(log(A(1:k)./A(k)));
    Ca(k) = k/length(A) * A(k)^Ha(k);
    Qa(k) = (1.0e-2/Ca(k))^(-1/Ha(k));
end

plot(1:length(Hg), Qg, 1:length(Hl), Ql, 1:length(Ha), Qa, 'LineWidth', ...
     2);
grid on

% plot(1:length(Hg), Hg, 1:length(Hl), Hl, 1:length(Ha), Ha, 'LineWidth', ...
%      2);
% grid on

% legend('gains', 'losses', 'absolute values');
xlabel('k: # upper order statistics');
% ylabel('Tail index \alpha');
ylabel('99\% quantile');

% alfa = [2.72, 3.00, 2.89];
% K = [289, 232, 287];
% N = [10389, 8871, 19260];
% Q = [gains(K(1))^alfa(1), losses(K(2))^alfa(2), A(K(3))^alfa(3)];
% C = K ./ N .* Q;

% q99 = [(1.0e-2/C(1))^alfa(1), (1.0e-2/C(2))^alfa(2),...
%        (1.0e-2/C(3))^alfa(3)];
