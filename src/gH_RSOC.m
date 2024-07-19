%GH_RSOC implements a membership and barrier function oracle for the rotated
% second-order cone
%   { x \in R_+^2 x R^n | x(1)*x(2) >= ||x(3:end)||^2 }
% by a transformation to the second-order cone.
% It requires no parameters.
% --------------------------------------------------------------------------
% USAGE of "gH_RSOC"
% [in, g, Hi, Li] = gH_RSOC(x)
% --------------------------------------------------------------------------
% INPUT
% x:            primal iterate
%
% OUTPUT
% in:	0 if x is not in the interior of the cone. 1 if x is in the
%       interior of the cone.
% g:	gradient of the barrier function at x
% Hi:	function representing the inverse Hessian action at x
% Li:   function representing the inverse Cholesky action or similar
%
% The last 3 output may be anything if in==0.
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% gH_SOCP: membership and barrier function oracle for second-order cones.
% -------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz
% --------------------------------------------------------------------------

function [in, g, Hi, Li] = gH_RSOC(x, ~)

    n     = length(x);
    p     = x(1);
    q     = x(2);
    xbar  = x(3:end);

    in = p > 0 && q > 0 && 2*p*q > xbar'*xbar;
    
    if in
        if nargout == 1
            return;
        end
        
        T = blkdiag([sqrt(1/2), sqrt(1/2); sqrt(1/2), -sqrt(1/2)], speye(n-2)); % note: equals its inverse
        switch nargout
            case 2
                [~, g] = gH_SOCP(T*x, []);
            case 3
                [~, g, Hi] = gH_SOCP(T*x, []);
            case 4
                [~, g, Hi, Li] = gH_SOCP(T*x, []);
        end
        
        g = T*g;
        
        if nargout > 2
            %H = T*H*T;
            Hi = @(v)(T*Hi(T*v));
        end
        
        if nargout > 3
            %L = T*L;
            Li = @(M)(Li(T*M));
        end
          
    else
        g = NaN; Hi = NaN; Li = NaN;
    end
    
return
