%gH_SOCP implements the usual logarithmic barrier function for the nonnegative orthant

function [in, g, H, L] = gH_LP(x, ~)

    n  = length(x);
    in = min(x)>0;
    
    if in
        g = -1./x;
        if nargout > 2
            H = sparse(1:n,1:n,x.^(-2),n,n,n);
            L = sparse(1:n,1:n,-g,n,n,n);
        end
    else
        g = NaN; H = NaN; L = NaN;
    end

return
