clear all
close all

fold = 40000;
% fold = 1;



% alpha = 2.17
% params = [
%     0.1068,  0.8923
    % 0.2009, 0.7962
    % 0.2991,  0.6946
    % 0.4001, 0.5889
    % 0.4996,  0.4835
    % 0.5995, 0.3786
    % 0.6992,  0.2717
    % 0.8009, 0.1608
    % 0.8970,  0.0541
    %        ];

% alpha = 2.40;
% params = [
%     0.0663, 0.9329 %, a = 2.3951
%     0.1017, 0.8964 % , a = 2.4019
%     0.2016, 0.7914 % , a = 2.3968
%     0.3013, 0.6838 % , a = 2.3955
%     0.4007, 0.5735 % , a = 2.3960
%     0.5011, 0.4593 % , a = 2.3951
%     0.5992, 0.3443 % , a = 2.3973
%     0.6992, 0.2221 % , a = 2.4045
%     0.8013, 0.0932 % , a = 2.4047    
%     0.8684, 0.0040 % , a = 2.4048
%          ];

% alpha = 3.00;
% params = [
%     0.1041, 0.8908 % a = 2.9902
%     0.1994, 0.7822 % , a = 3.0013
%     0.3008, 0.6578 %, a = 3.0097
%     0.4000, 0.5273 %, a = 3.0065
%     0.5000, 0.3874 % a = 2.9936, tau = 0.9999
%     0.6031, 0.2259 % a = 3.0088, tau = 1.9996    
%     0.7156, 0.0289 % a = 3.0099, tau = 2.2151
%          ];

% alpha = 4.00;
% params = [
%     0.1012, 0.8885 % , a = 4.0097
%     0.1998, 0.7592 % a = 4.0080
%     0.2996, 0.6060 % , a = 4.0027, tau = -4.1641    
%     0.4001, 0.4239 %, a = 4.0037
%     0.5004, 0.2082 % , a = 3.9908, tau = 1.0023
%     0.5491, 0.0784 % a = 4.0086, tau = 1.2009    
%          ];

% alpha = 5.00;
% params = [
%     0.1042, 0.8784 % , a = 4.9994
%     0.2008, 0.7302 %, a = 4.9974
%     0.2999, 0.5332 %, a = 4.9958
%     0.3999, 0.2681 %, a = 5.0035    
%     0.4497, 0.1009 %, a = 5.0042
%          ];

% Fixed alpha
% params = [
%     0.1000, 0.8300 % , a = 10.6775    
%     0.1000,  0.840 % 0, a = 9.9250    
%     0.1000,  0.850 % 0, a = 9.0759    
%     0.1000,  0.860 % 0, a = 8.0930    
%     0.1000,  0.870 % 0, a = 6.9546    
%     0.1000,  0.880 % 0, a = 5.6097    
    
%     0.1500,  0.780 % 0, a = 6.7646
%     0.1500,  0.790 % 0, a = 6.2856    
%     0.1500,  0.800 % 0, a = 5.7572    
%     0.1500,  0.810 % 0, a = 5.1741    
%     0.1500,  0.820 % 0, a = 4.5221    
%     0.1500,  0.830 % 0,  = 3.7885    

%     0.2000,  0.730 % 0, a = 5.0468    
%     0.2000,  0.740 % 0, a = 4.7146
%     0.2000,  0.750 % 0, a = 4.3551    
%     0.2000,  0.760 % 0, a = 3.9663    
%     0.2000,  0.770 % 0, a = 3.5420    
%     0.2000,  0.780 % 0, a = 3.0775    
    
%     0.2500,  0.680 % 0, a = 4.1305
%     0.2500,  0.690 % 0, a = 3.8860    
%     0.2500,  0.700 % 0, a = 3.6252    
%     0.2500,  0.710 % 0, a = 3.3462    
%     0.2500,  0.720 % 0, a = 3.0470    
%     0.2500,  0.730 % 0, a = 2.7249    
    
%     0.3000,  0.630 % 0, a = 3.5802
%     0.3000,  0.640 % 0, a = 3.3924    
%     0.3000,  0.650 % 0, a = 3.1940    
%     0.3000,  0.660 % 0, a = 2.9838    
%     0.3000,  0.670 % 0, a = 2.7608    
%     0.3000,  0.680 % 0, a = 2.5235    
%          ];


% dist = struct('name', 'Garch1_1', 'prmt', [2.3e-6, 0.15, 0.84], 'distr', ...
%               struct('Name', 'Gaussian'), 'TailExponent', 2.9664);

params =  [
    0.1068, 0.8923, 2.17
    0.1017, 0.8964, 2.40
    0.1041, 0.8908, 3.00
    0.1012, 0.8885, 4.00
          ];

for q = [24, 10, 50]
    for N = [10]
        for n = 1 : size(params, 1)
            a1 = params(n, 1);
            b1 = params(n, 2);
            a = params(n, 3);
            if (exist(sprintf(['../data/tail_exponent_%.2f/'...
                              'GarchWishartOrdN%dQ%dA%.4fB%.4fEig.mat'], ...
                              a, N, q, a1, b1), 'file') == 2)
                
                load(sprintf(['../data/tail_exponent_%.2f/'...
                              'GarchWishartOrdN%dQ%dA%.4fB%.4fEig.mat'], ...
                             a, N, q, a1, b1), 'ev');
                load(sprintf(['../data/tail_exponent_%.2f/' ...
                              'GarchWishartOrdN%dQ%dA%.4fB%.4fX.mat'], ...
                             a, N, q, a1, b1), 'X');
                load(sprintf(['../data/tail_exponent_%.2f/' ...
                              'GarchWishartOrdN%dQ%dA%.4fB%.4fC.mat'], ...
                             a, N, q, a1, b1), 'C');
            end
            dist = struct('name', 'Garch1_1', 'prmt', [1- a1 - b1, a1, b1], 'distr', ...
                          struct('Name', 'Gaussian'), 'TailExponent', ...
                          2.30);
            T = N * q;
            ev1 = NaN(N, fold);
            X1 = NaN(N, N, fold);
            C1 = NaN(N, N, fold);
            for i = 1:fold
                R = gen_ret_mtx(N, T, dist, 0, []);
                R = diag(1./std(R')) * R;
                C1(:, :, i) = R*R' ./ T;
                [V, D] = eig(C1(:, :, i));
                ev1(:, i) = D(logical(eye(N)));
                X1(:, :, i) = V;
            end
            % ev = ev1;
            % X = X1;
            % C = C1;
            ev = [ev, ev1];
            C(:, :, end+1 : end+fold) = C1;
            X(:, :, end+1 : end+fold) = X1;

            save(sprintf(['../data/tail_exponent_%.2f/'...
                         'GarchWishartOrdN%dQ%dA%.4fB%.4fEig.mat'], ...
                         a, N, q, a1, b1), 'ev');
            save(sprintf(['../data/tail_exponent_%.2f/' ...
                         'GarchWishartOrdN%dQ%dA%.4fB%.4fX.mat'], ...
                         a, N, q, a1, b1), 'X');
            save(sprintf(['../data/tail_exponent_%.2f/' ...
                         'GarchWishartOrdN%dQ%dA%.4fB%.4fC.mat'], ...
                         a, N, q, a1, b1), 'C');
            fprintf('Saved files %.4f, %.4f.\n', a1, b1);
            
        end
    end
end
quit