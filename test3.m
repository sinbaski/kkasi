clear all
for n = 0:19
    intra_autocorr('Nordea_Bank', '2012-04-20', '2012-04-20', ...
                   30, 0, 0);
    fprintf('n=%d\n', n);
end
