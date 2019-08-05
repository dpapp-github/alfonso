%gH_RSOC implements a barrier function for the rotated second-order cone
%using whatever barrier for the (not rotated) second order cone is available

function [in, g, H, L] = gH_RSOC(x, ~)

    n     = length(x);
    p     = x(1);
    q     = x(2);
    xbar  = x(3:end);

    in = p > 0 && q > 0 && 2*p*q > xbar'*xbar;
    
    if in
        if nargout == 1
            return;
        end
        
        T = blkdiag([sqrt(1/2), sqrt(1/2); sqrt(1/2), -sqrt(1/2)], speye(n-2));
        
        switch nargout
            case 2
                [~, g] = gH_SOCP(T*x, []);
            case 3
                [~, g, H] = gH_SOCP(T*x, []);
            case 4
                [~, g, H, L] = gH_SOCP(T*x, []);
        end
        
        g = T*g;
        
        if nargout > 2
            H = T*H*T;
        end
        
        if nargout > 3
            L = T*L;
        end
          
    else
        g = NaN; H = NaN; L = NaN;
    end
    
return
