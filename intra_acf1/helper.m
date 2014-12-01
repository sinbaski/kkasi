clear all

% intra_autocorr('Nordea_Bank', '2012-01-16', '2012-03-15', 30, 1);
% intra_autocorr('Nordea_Bank', '2012-01-16', '2012-03-15', 30, 0, 0);

% % intra_autocorr('Nordea_Bank', '2012-01-16', '2012-03-15', 20, 1);
% intra_autocorr('Nordea_Bank', '2012-01-16', '2012-03-15', 20, 0);

% intra_autocorr('Nordea_Bank', '2012-01-16', '2012-03-15', 10, 0);
% %intra_autocorr('Nordea_Bank', '2012-01-16', '2012-03-15', 10, 0);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intra_autocorr('Nordea_Bank', '2012-03-16', '2012-04-20', 30, 1, 0);
% % intra_autocorr('Nordea_Bank', '2012-03-16', '2012-04-20', 30, 0);

% intra_autocorr('Nordea_Bank', '2012-03-16', '2012-04-20', 20, 0, 1);
% % intra_autocorr('Nordea_Bank', '2012-03-16', '2012-04-20', 20, 0);

% intra_autocorr('Nordea_Bank', '2012-03-16', '2012-04-20', 10, 0, 1);
% intra_autocorr('Nordea_Bank', '2012-03-16', '2012-04-20', 10, 0);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% intra_autocorr('Ericsson_B', '2012-01-25', '2012-02-28', 30, 1);
% intra_autocorr('Ericsson_B', '2012-01-25', '2012-02-28', 30, 0);

% intra_autocorr('Ericsson_B', '2012-01-25', '2012-02-28', 20, 1);
% intra_autocorr('Ericsson_B', '2012-01-25', '2012-02-28', 20, 0);

% intra_autocorr('Ericsson_B', '2012-01-25', '2012-02-28', 10, 1);
% intra_autocorr('Ericsson_B', '2012-01-25', '2012-02-28', 10, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intra_autocorr('Nordea_Bank', '2013-09-03', '2013-09-03', ...
                   5, 0, 0);
for n = 1:40
intra_autocorr('Nordea_Bank', '2013-09-03', '2013-09-03', ...
                   5, n, 0);
    fprintf('Finished round %d\n', n);
end

