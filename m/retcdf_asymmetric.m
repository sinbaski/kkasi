function p = retcdf_asymmetric(x, psi, sigma, evbar)
% N = 5;
% x = x./evbar;
% Sk = zeros(size(x));
% for k =0:N
%     A = (-1)^k/(2*k+1)/factorial(k)/(2*(1-psi^2))^(k+1/2);
%     Si = zeros(size(x));
%     for i = 0:2*k+1
%         B = nchoosek(2*k+1, i)*psi^(2*k+1-i).*x.^i.*...
%             exp(sigma^2*i^2/2);
%         Sj = zeros(size(x));
%         for j = 0:2*k+1-i
%             C = nchoosek(2*k+1-i, j)*2^(j/2)*(sigma*i)^(2*k+1-i-j)*...
%                 sqrt(2)*gamma((j+1)/2);
%             Sj = Sj + C;
%         end
%         B = B * Sj;
%         Si = Si + B;
%     end
%     A = A * Si;
%     Sk = Sk + A;
% end
% p = 1/2 - 1/sqrt(2)/pi*Sk;

x = x./evbar;
u1 = 0.27;
u2 = -1.1;
a1 = fzero(@(a) a - e^(-a*sigma)*a)