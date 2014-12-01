classdef nordea_15min_profile < handle
    properties (Access = public)
        s = 33
        lv = []
        w = []
        y = []
        MA = []
        % The Johnson Su structure with name-value pairs.
        jsp = []
    end
    
    methods (Access = public)
        % pdf of volatility at horizon h
        function p = sigma_pdf(me, v, h)
        %% conditional probability density of the volatility
        % v: vector or scalar The volatility values for which the
        % probability density is to be evaluated.
            p = zeros(size(v));
            v = v(v > 0);
            if isempty(v)
                return;
            end
            s = me.s;
            t = length(me.lv) + h;
            tw = t - s - 1;
            yh = log(v) - me.lv(t-1) - me.lv(t-s)...
                 + me.lv(t-s-1) + me.MA ...
                 * expectation(me.y, tw-[1:s+1]);
            p = johnson_su_pdf(me.jsp, yh);
            p = p ./ v;
        end
        
        function p = ret_pdf(me, r, h)
            p = NaN(size(r));
            for n = 1:length(r)
                func = @(v) normpdf(r(n)./v, 0, 1) .* ...
                       me.sigma_pdf(v, h) ./ v;
                p(n) = integral(func, 0, Inf);
            end
        end
        
        function c = ret_cdf(me, r, h)
            c = NaN(size(r));
            for n = 1:length(r)
                if r(n) == -Inf
                    c(n) = 0;
                else
                    c(n) = integral(@(x) me.ret_pdf(x, h), -Inf, r(n));
                end
            end
        end
        
        function [w_f, lv_f] = forecast(me, h)
            n = length(me.w);
            m = length(me.lv);
            s = me.s;
            for j = 1:h
                me.w(n+j) = -me.MA * expectation(me.y, n+j-[1:s+1]);
                me.lv(m+j) = me.w(n+j) + me.lv(m+j-1) ...
                    + me.lv(m+j-s) - me.lv(m+j-s-1);
            end
            w_f = me.w(n+1 : end);
            lv_f = me.lv(m+1 : end);
            me.w = me.w(1:n);
            me.lv = me.lv(1:m);
        end
        
        function obj = update(me, lv_new)
            m = length(me.lv);
            n = length(me.w);
            l = length(me.y);
            s = me.s;
            me.lv = [me.lv; lv_new];
            for k = 1:length(lv_new)
                me.w(n+k) = me.lv(m+k) - me.lv(m+k-1) - me.lv(m+k-s) ...
                    + me.lv(m+k-s-1);
                me.y(l+k) = me.w(n+k) + me.MA*me.y(l+k-[1:s+1]);
            end
            obj = me;
        end
        
        function res = simulate_residual(me, N, np)
            z = randn(N+(s+1), np);
            res = me.jsp.Xi + me.jsp.Lambda * ...
                sinh((z - me.jsp.Gamma)./me.jsp.Delta);
        end
        
        function obj = nordea_15min_profile(lv, w, y, MA, jsp)
            obj.lv = lv;
            obj.w = w;
            obj.y = y;
            obj.MA = MA;
            obj.jsp = jsp;
        end
    
    end
end
