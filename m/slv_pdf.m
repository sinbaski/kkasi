function y = slv_pdf(x, psi, sig)
x = reshape(x , length(x), 1);
da = 1.0e-2;
a = -10:da:10;
Ki = - 2 * psi .* x * (a .*  exp(-a .* sig))...
     + x.^2 * exp(-2 * a * sig);
Ki = -Ki ./ 2 ./ (1- psi^2);
K = exp(Ki);

Kr = exp(-((a + sig * (1 - psi^2)).^2 - sig^2 * (1 - psi^2)^2) ./...
         (2 * (1 - psi^2))) .* da;

y = K * Kr' ./ 2 ./ pi ./ (1 - psi^2);

