%% log-likelihood function of the model
% w_t = (1 - theta B)(1 - Theta B^s) alpha_t
% s fixed to 33;

function l = loglikelihood(param, model, MALags, obs)
n = length(obs);
s = 33;

% a: gamma
% b: delta
% c: lambda
% m: xi
% a = param(1);
% b = param(2);
% c = param(3);
% m = param(4);
% theta = param(1);
% Theta = param(2);
% mmt = johnson_su_moments12(param);

% y = ma_infer(obs, theta, Theta, s);
% y = ma_infer(obs, param, MALags, s);
c = 1;
for k = MALags
    model.MA{k} = param(c);
    c = c + 1;
end
y = infer(model, obs);
moments = [mean(y), var(y), skewness(y), kurtosis(y)];
if (moments(3) < 0 || moments(3) > 0.5 ||...
    moments(4) < 3.6 || moments(4) > 8.0)
    l = -Inf;
    return;
end

jsp = johnson_su_params(moments);
a = jsp(1);
b = jsp(2);
c = jsp(3);
m = jsp(4);

z = b*asinh((y - m)/c) + a;
l = n*log(b/c) - z'*z/2  - sum(log(1 + ((y - m)/c).^2))/2 - n* ...
    log(2*pi)/2;

% gradz = johnson_su_grad(param([1:4]), y, z);
% grad = NaN(4, 1);
% grad(1) = 0;
% grad(2) = -z' * gradz(:, 2) + n/b;
% grad(3) = sum([c^3 ./ (y - m).^2 + c].^(-1)) - z' * gradz(:, 3) - n/c;
% grad(4) = sum((y - m) ./ (c^2 + (y - m).^2)) - z' * gradz(:, 4);
