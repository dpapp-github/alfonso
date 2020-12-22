%GH_RSOC implements a membership and barrier function oracle for rotated
% second-order cones, by transformation to the second-order cone.
% It requires no parameters.
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
% USAGE of "gH_RSOC"
% [in, g, H, L] = gH_RSOC(x)
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
% gH_SOCP: membership and barrier function oracle for second-order cones.
% -------------------------------------------------------------------------

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
