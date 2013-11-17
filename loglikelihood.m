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
theta1 = param(5);
theta2 = param(6);

coef = zeros(1, s+1);
coef([1, s, s+1]) = [theta1, theta2, theta1*theta2];

y = ma_infer(obs, coef);
z = b*asinh((y - m)/c) + a;

l = -0.5*z'*z -0.5*sum(log(1 - ((y - m)/c).^2));


