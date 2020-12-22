%GH_SOCP implements a membership and barrier function oracle for second-
% order cones, using the standard logarithmic barrier. It requires no
% parameters.
% --------------------------------------------------------------------------
%
% Copyright (C) 2018-2020 David Papp and Sercan Yildiz.
%
% Redistribution and use of this software are subject to the terms of the
% 2-Clause BSD License. You should have received a copy of the license along
% with this program. If not, see <https://opensource.org/licenses/BSD-2-Clause>.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz    <syildiz@qontigo.com>  
%
% --------------------------------------------------------------------------
% USAGE of "gH_SOCP"
% [in, g, H, L] = gH_SOCP(x)
% --------------------------------------------------------------------------
% INPUT
% x:            primal iterate
%
% OUTPUT
% in:	0 if x is not in the interior of the cone. 1 if x is in the
%       interior of the cone.
% g:	gradient of the barrier function at x
% H:	Hessian of the barrier function at x
% L:	Cholesky factor of the barrier function at x
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

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
