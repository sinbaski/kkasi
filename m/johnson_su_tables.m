function [Gamma, Delta] = johnson_su_tables(m3, m4)
sk = 0:0.05:0.5;
kts = [3.2:0.1:8.0, 8.2:0.2:9.0];

if m3 < min(sk) || m3 > max(sk) || m4 < min(kts) || m4 > max(kts)
    Gamma = NaN;
    Delta = NaN;
    return;
end
        
[v, x] = min(abs(sk - m3));
[v, y] = min(abs(kts - m4));
gammas = importdata('johnson_su_gamma.txt');
deltas = importdata('johnson_su_delta.txt');
Gamma = -gammas(y, x);
Delta = deltas(y, x);
