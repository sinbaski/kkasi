%% log-likelihood function of the model
% w_t = (1 - theta B)(1 - Theta B^s) alpha_t
% s fixed to 33;

function l = loglikelihood(param, obs)
n = length(obs);
s = 33;

a = param(1);
b = param(2);
c = param(3);
m = param(4);
theta = param(5);
Theta = param(6);
% mmt = johnson_su_moments12(param);

y = ma_infer(obs, theta, Theta, s);

z = b*asinh((y - m)/c) + a;
l = -z'*z/2 -sum(log(1 + ((y - m)/c).^2))/2;

% gradz = johnson_su_grad(param([1:4]), y, z);
% grad = NaN(4, 1);
% grad(1) = 0;
% grad(2) = -z' * gradz(:, 2) + n/b;
% grad(3) = sum([c^3 ./ (y - m).^2 + c].^(-1)) - z' * gradz(:, 3) - n/c;
% grad(4) = sum((y - m) ./ (c^2 + (y - m).^2)) - z' * gradz(:, 4);
