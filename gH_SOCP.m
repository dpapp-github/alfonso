%gH_SOCP implements the usual logarithmic barrier function for the second-order cone

function [in, g, H, L] = gH_SOCP(x, ~)

    n     = length(x);
    xbar  = x(2:end);
    det   = x(1)^2 - xbar'*xbar;

    in = x(1) > 0 && det > 0;
    
    if in
        g = 2*[-x(1); xbar] / det;
        
        if nargout > 2
            H = diag( 2*[-1/det; repmat(1/det,n-1,1)] ) + g*g';
        end
        
        if nargout > 3            
            R = diag( sqrt(2*repmat(1/det,n,1)) );
            R = cholupdate(R, g);
            [R, err] = cholupdate(R, [2/sqrt(det); zeros(n-1,1)], '-');
            L = R';

            if err > 0
                in = false; g = NaN; H = NaN; L = NaN;
                return;
            end
        end
        
    else
        g = NaN; H = NaN; L = NaN;
    end
    
return
