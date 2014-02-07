clear all
syms N A b;

assume(N, 'integer');
assume(b > 0);
assumeAlso(b < 1);

R = sym('r%d%d', [2, 2]);
assume(R, 'real');

M = eye(3);
% b = 0.5;
M = sym(M);
for i = 1:size(M, 1)-1
    M(i, i+1) = -b;
end

% A = R*M;

% [V1, E1] = eig(R*R');
% [V2, E2] = eig(R'*R);
[V, E] = eig(M*M');

% t1 = trace(A*A');
% t2 = trace(A*A'*R*R');

% isequaln(R'*V1, V2)
