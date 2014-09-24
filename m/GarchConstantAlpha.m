clear all
f = @(a, b) a*b;
delf1 =@(a, b, c) a;
delf1 =@(a, b, c) b;
delf = @(a, b, c) [delf1(a, b, c); delf2(a, b, c)];

alpha = 3;

n = 4;
a1 = 0.55;
b1 = 0.08;

% Matrix for rotating pi/2 clockwise
RC = [ 0,  1
      -1,  0];

RCC = [0,  -1
       1,  0];

V = [a1; b1];
% Search along larger b1 and smaller a1
% Rotate the gradient of alpha clockwise.
p = alpha;
q = NaN;
for k = 1:n
    % Calculate the gradient of alpha
    G = delf(a1, b1, p);
    PG = RC * G;
    PG = PG ./ norm(PG);
    V = V + PG .* norm(V)/20;
    
    % Calculate the new alpha
    q = f(V(1), V(2));
    if abs(q - alpha)/alpha < 0.001
        a1 = V(1);
        b1 = V(2);
        p = q;
        fprintf('a1 = %.3f, b1 = %.3f, alpha = %.3f\n', a1, b1, q);
        continue;
    end
    % if the new alpha has deviated too much for alpha,
    % Move along the gradient or its opposite.
    G = delf(V(1), V(2), q);
    Gnorm = norm(G);
    g = G ./ Gnorm;
    delta = abs(q - alpha)/Gnorm;

    if (q < alpha)
        V = V + delta * g;
    else
        V = V - delta * g;
    end
    a1 = V(1);
    b1 = V(2);
    q = f(V(1), V(2));
    fprintf('a1 = %.3f, b1 = %.3f, alpha = %.3f\n', a1, b1, q);
    p = q;
end
