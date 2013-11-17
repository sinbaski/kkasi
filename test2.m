clear all
% z = [
%     460
%     457
%     452
%     459
%     462
%     459
%     463
%     479
%     493
%     490
%     ];
% theta = 0.5;
% w = z(2:end) - z(1:end-1);
% a = ma_infer(w, theta);
func = @(x, y, z) x.*y.*z.^34 - x.*z - y.*z.^33 + 1;

