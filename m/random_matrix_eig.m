%% Generate random matrices
%function lambda = random_matrix_eig(N, T, Sigma, filename)
function lambda = random_matrix_eig(T, Sigma, filename)
% M = random('Normal', 0, 1, [N, T]);
% C = M * M'/T;
C = wishrnd(Sigma, T);
lambda = eig(C);
if exist(filename, 'file')
    save(filename, 'lambda', '-append', '-ascii', '-double', '-tabs');
else
    save(filename, 'lambda', '-ascii', '-double', '-tabs');
end
