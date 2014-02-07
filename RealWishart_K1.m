function [nf, df] = RealWishart_K1(n1, n2)
if mod(n1, 2) == 0 && mod(n2, 2) == 1
    nf = []; % factors of the numerator
    df = []; % factors of the denominator
    for i = 1:n1
        nf = [nf, 1:(n2-n1-1)/2+i-1];
        if mod(i, 2) == 0
            df = [df, 1:(n2-(i-1))/2-1, 1/2:(n1-(i-1))/2-1];
        else
            df = [df, 1/2:(n2-(i-1))/2-1, 1:(n1-(i-1))/2-1];
        end
    end
    u = length(nf) - length(df);
    if u > 1
        df = [ones(1, u), df];
    elseif u < -1
        nf = [ones(1, -u), nf];
    end
    % K1 = prod(nf ./ df);
elseif mod(n1, 2) == 0 && mod(n2, 2) == 0
    nf = []; % factors of the numerator
    df = []; % factors of the denominator
    pis = 0;
    for i = 1:n1
        nf = [nf, sqrt(pi), 1/2:(n2-n1-1)/2+i-1];
        if mod(i, 2) == 0
            df = [df, 1/2:(n2-(i-1))/2-1, 1/2:(n1-(i-1))/2-1];
        else
            df = [df, 1:(n2-(i-1))/2-1, 1:(n1-(i-1))/2-1];
        end
    end
    u = length(nf) - length(df);
    if u > 1
        df = [ones(1, u), df];
    elseif u < -1
        nf = [ones(1, -u), nf];
    end
    % K1 = prod(nf ./ df);
elseif mod(n1, 2) == 1 && mod(n2, 2) == 1
    nf = []; % factors of the numerator
    np = n1/2;
    df = []; % factors of the denominator
    for i = 1:n1
        nf = [nf, sqrt(pi), 1/2:(n2-n1-1)/2+i-1];
        if mod(i, 2) == 1
            df = [df, 1/2:(n2-(i-1))/2-1, 1/2:(n1-(i-1))/2-1];
        else
            df = [df, 1:(n2-(i-1))/2-1, 1:(n1-(i-1))/2-1];
        end
    end
    nf = [nf, 1/2:(n2+n1+1)/2-1];
    nf = [ones(1, (n1+n2)/2+1)*2, nf];
    df = [sqrt(2), df];

    u = length(nf) - length(df);
    if u > 1
        df = [ones(1, u), df];
    elseif u < -1
        nf = [ones(1, -u), nf];
    end
    % K1 = prod(nf ./ df);
elseif mod(n1, 2) == 1 && mod(n2, 2) == 0
    nf = []; % factors of the numerator
    df = []; % factors of the denominator
    for i = 1:n1
        nf = [nf, 1:(n2-n1-1)/2+i-1];
        if mod(i, 2) == 0
            df = [df, 1/2:(n2-(i-1))/2-1, 1:(n1-(i-1))/2-1];
        else
            df = [df, 1:(n2-(i-1))/2-1, 1/2:(n1-(i-1))/2-1];
        end
    end
    nf = [nf, 1:(n2+n1+1)/2-1];
    nf = [ones(1, (n1+n2+1)/2)*2, nf];

    u = length(nf) - length(df);
    if u > 1
        df = [ones(1, u), df];
    elseif u < -1
        nf = [ones(1, -u), nf];
    end
    % K1 = prod(nf ./ df);
end
