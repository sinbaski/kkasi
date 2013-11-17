%% Infer innovations of a MA process
function a = ma_infer(w, theta)
n = length(w);
q = length(theta);

a = zeros(3*q+n, 1);
a_off = 2*q;

e = zeros(n+3*q, 1);
e_off = q;

w = [zeros(q, 1); w; zeros(q, 1)];
w_off = q;

b = NaN(n, 1);
% number of iterations already performed
it = 0;

while it == 0 || mean(abs(a(a_off + [1:n]) - b)) > 1e-4;
    % record the current a for later comparison
    b = a(a_off + [1:n]);
    
    % update e from n+q to 1, using w from n to n+q
    % e from n+q+1 to n+2q are assumed 0.
    % e[t] = w[t] + sum(i=1:q) e[t+i];
    for t = n+q : -1 : 1
        e(t + e_off) = w(t + w_off) + theta * e(t + e_off + [1:q]);
    end

    % back forecast w from 0 to -q+1
    % E(e[t]) = E(w[t]) + sum{j=1:q} theta[j]*E(e[t+j]);
    % E(e[t]) = 0 for t <= 0 hence
    % E(w[t]) = -sum{j=1:q} theta[j]*E(e[t+j])
    for t = 0:-1:1-q
        w(t + w_off) = - theta * e(t + e_off + [1:q]);
    end

    % forecast a from 1-q to n
    % E(a[t]) = E(w[t]) + sum{i=1:q} E(a[t-i]) * theta[i]
    % a from -2q+1 to -q are assumed 0
    for t = 1-q : n
        a(t + a_off) = w(t + w_off) + theta * a(t + a_off - [1:q]);
    end

    % forecast w from n+1 to n+q using infered
    % innovations from the last iteration
    for t = n+1 : n+q
        w(t + w_off) = a(t + a_off) - theta * a(t + a_off - [1:q]);
    end

    it = it + 1;
end
