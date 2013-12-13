clear all
close all
load('/tmp/models.mat');

suppress = @(phi, t, tau) exp(-t./tau) .* phi;

T = 1e3;
N = 50;
fold = 4e3;
n = (N^2 - N)/2;
%Cij = NaN(1, n * fold);
%u = 0;
colors = [
    0   0     1
    0   1/8   7/8,
    0   2/8   6/8
    0   3/8   5/8
    0   4/8   4/8
    0   5/8   3/8
    0   6/8   1/8
    0   1     0
    ];

tau = [1:5, [1:6]*4, Inf];
gen_data = zeros(1, length(tau));
for l = 9:length(tau)
    if (exist(sprintf('./eigen-%d.mat', l), 'file') ~= 2 ||...
        gen_data(l) == 1)
        
        ev_new = NaN(1, N*fold);
        model = model0;
        model.AR = num2cell(suppress(...
            cell2mat(model0.AR),...
            [1:length(model.AR)],...
            tau(l)));
        model.SAR = num2cell(suppress(...
            cell2mat(model0.SAR),...
            [1:length(model.SAR)],...
            tau(l)));
        for k = 1:fold
            [LV, E, V] = simulate(model, T, 'numPaths', N);
            M = exp(LV) .* randn(size(LV));
            sigmas = diag(1./std(M));
            % Normalize the simulated returns
            M = M*sigmas;
            C = M'*M/T;
            ev_new((k-1)*N+1 : k*N) = eig(C);
            % for l=1:N
            %     for m = l+1:N
            %         Cij(u + 1) = C(l, m);
            %         u = u + 1;
            %     end
            % end
        end
        if exist(sprintf('./eigen-%d.mat', l), 'file') == 2
            load(sprintf('./eigen-%d.mat', l));
            ev = [ev, ev_new];
        else
            ev = ev_new;
        end
        save(sprintf('./eigen-%d.mat', l), 'ev');
    else
        load(sprintf('./eigen-%d.mat', l));
    end
    % x = linspace(min(ev), max(ev), 40);
    % count = hist(ev, x);
    % plot(x, count/length(ev)/(x(2)-x(1)), 'Color', colors(l, :));
    %% Investigate the distribution of the largest eigen value
    
end
% I1 = ev < 0.47;
% I2 = ev > 1.6;

% q = N/T;
% x = [(1 - sqrt(q))^2 : 1e-3: (1 + sqrt(q))^2];
% % x = [2.7 : 1e-3: (1 + sqrt(q))^2];
% rho = MarcenkoPasturPDF([q, 1], x);
% plot(x, rho, 'r');
% hold off

