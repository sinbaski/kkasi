clear all
close all
timefmt = 'yyyy-mm-dd HH:MM:SS';
% in units of minute
interval = 5;
start_day = '2012-03-16';
end_day = '2012-04-20';
company = 'nordea_bank';

get_data = 0;

ret = get_intra_ret_simple(company, '2012-01-16', '2012-04-20', interval);
r = (ret - mean(ret))./std(ret);

q1 = [-2.1, -0.7, 0.7, 2.1];
probs = normcdf(q1);
q2 = quantile(r, probs);
[R1, type1] = johnsrnd([q1; q2], 1e5, 1);
[F1, X1] = ecdf(R1);

ecdf(r);
hold on
stairs(X1, F1, 'r');
%stairs(X2, F2, 'g');
grid on
hold off;


