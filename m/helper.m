clear all
% alpha = 4.00
params = [
    0.1012, 0.8885 % , a = 4.0097
    0.2996, 0.6060 % , a = 4.0027, tau = -4.1641    
    0.5004, 0.2082 % , a = 3.9908, tau = 1.0023
    0.5491, 0.0784 % a = 4.0086, tau = 1.2009    
         ];

for n = 1 : size(params, 1)
    fprintf(['mv GarchWishartOrdN50Q24A%.4fB%.4fEig.mat ' ...
             'tail_exponent_4.00/GarchWishartOrdN50Q24A%.4fB%.4fEig.mat\n'],...
            params(n, 1), params(n, 2), params(n, 1), params(n, 2));
    fprintf(['mv GarchWishartOrdN50Q24A%.4fB%.4fC.mat ' ...
             'tail_exponent_4.00/GarchWishartOrdN50Q24A%.4fB%.4fC.mat\n'],...
            params(n, 1), params(n, 2), params(n, 1), params(n, 2));
    fprintf(['mv GarchWishartOrdN50Q24A%.4fB%.4fX.mat ' ...
             'tail_exponent_4.00/GarchWishartOrdN50Q24A%.4fB%.4fX.mat\n'],...
            params(n, 1), params(n, 2), params(n, 1), params(n, 2)); 
end
