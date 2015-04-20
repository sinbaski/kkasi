% The Green's function
% Z: Column vector
% G: Column vector
function G = LognormalGreen(Z, v, q, varargin)
    a = (1 - sqrt(q))^2;
    b = (1 + sqrt(q))^2;
    MPGreen = @(z) (z + q - 1 - sqrt((z-a).*(z-b)))./(2*q.*z);
    % g = MPGreen(Z);
    G = NaN(length(Z), 1);
    options = optimset('Jacobian','on', 'Display', 'off');
    if (length(varargin) > 0)
        Ginit = varargin{1};
    else
        Ginit = [];
    end
    for n = 1 : length(Z)
        % if (Z(n) > a && Z(n) < b)
        %     initial = MPGreen(Z(n));
        % else
        %     initial = 0;
        % end
        if ~isempty(Ginit)
            initial = Ginit(n);
        elseif n > 1
            initial = G(n-1);
        elseif Z(n) > a && Z(n) < b
            initial = MPGreen(Z(n));
        else
            initial = 0;
            %        else
            % g = MPGreen(Z(n));
            % if abs(imag(G(n-1))) < abs(imag(g))
            %     initial = g;
            % else
            %     initial = G(n-1);
            % end
        end
        [U, val] = lsqnonlin(@(z) LognormalGreenFun1(z(1) + i*z(2),...
                                                     v, q, [Z(n); ...
                            0]), [real(initial); imag(initial)], [-Inf, ...
                            -Inf], [Inf, -1.0e-6], options);
        if U(2) >= 0
            break;
        end
        G(n) = U(1) + i*U(2);
    end
    if n < length(Z)
        G(n:end) = NaN;
    end
