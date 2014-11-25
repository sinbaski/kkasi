% spectral density of the wishart matrix
clear all

n = 500; % number of time series
T = 2000;
distr.name = 'normal';
generate_data = 0;

tau = 1;
% distr.name = 'lorentzian';
% distr.name = 'student t';
% distr.prmt = 4;

%auto_corr = @(t) exp(-t/tau);
auto_corr = @(t) 0;

filename = 'data/wishart_eig.txt';

if generate_data
    for l = 1:2000
        R = gen_ret_mtx(n, T, distr, auto_corr);
        %C = cpt_cov(R, distr);
        C = (1/T) * R * R';
        lambda = eig(C);
        if exist(filename, 'file')
            save(filename, 'lambda', '-append', '-ascii', '-double', '-tabs');
        else
            save(filename, 'lambda', '-ascii', '-double', '-tabs');
        end
    end
else
    lambda = load(filename);

    cen = linspace(min(lambda), max(lambda), 400);
    dx = cen(2) - cen(1);
    empirical = double(hist(lambda, cen)) ./ (length(lambda) * dx);

    % [q, sigma]
    %prmt0 = [0.25, 1];
    b = max(lambda) * 1.1;
    a = min(lambda) * 0.9;


    q0 = (sqrt(b) - sqrt(a))^2/(sqrt(a) + sqrt(b))^2;
    var0 = (a + b) / (2 + 2*q0);

    %phat = mle(lambda, 'pdf', @MarcenkoPasturPDF, 'start', [q0,
    %var0]);
    [prmt, resnorm] = lsqcurvefit(@MarcenkoPasturPDF, [q0, var0], cen, empirical);

    theoretical = MarcenkoPasturPDF(prmt, cen);
    plot(cen, empirical, 'bo', cen, theoretical, 'r-');
end
