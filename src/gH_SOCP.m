%GH_SOCP implements a membership and barrier function oracle for the
% second-order cone 
%    { x \in R^n | x(1) >= norm(x(2:end),2) },
% using the standard logarithmic barrier.
% --------------------------------------------------------------------------
% USAGE of "gH_SOCP"
% [in, g, Hi, Li] = gH_SOCP(x)
% --------------------------------------------------------------------------
% INPUT
% x:        primal iterate
%
% OUTPUT
% in:       0 if x is not in the interior of the cone. 1 if x is in the
%           interior of the cone.
% g:        gradient of the barrier function at x
% Hi:	    function representing the inverse Hessian action at x
% Li:       function representing the inverse Cholesky action or similar
%
% The last 3 output may be anything if in==0.
% --------------------------------------------------------------------------
% Details on the last two outputs: denoting the Hessian at x by H,
%  - Li(M) returns L\M, where L is any matrix satisfying LL' = inv(H)
% Needs to work for a matrix.
%  - Hi(v) returns H\v.
% Here, v can be assumed to be a column vector (not a matrix).
%
% It is often sensible, but not always necessary or most efficient, to
% implement Hi using Li or an explicitly computed Cholesky factor.
% --------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz
% --------------------------------------------------------------------------

function [in, g, Hi, Li] = gH_SOCP(x, ~)

    n     = size(x,1);
    xbar  = x(2:n);
    det   = x(1)^2 - xbar'*xbar;

    g = NaN; Hi = NaN; Li = NaN;  % will be overwritten if in==1

    in = x(1) > 0 && det > 0;
    
    if in
        g = 2*[-x(1); xbar] / det;
        
        if nargout >= 3
            Hi = @(v)( x*(x'*v) - det/2 * [v(1,:); -v(2:end,:)] );
            if nargout == 4
                R = diag( repmat(sqrt(2/det),n,1) );
                R = cholupdate(R, g);
                [R, err] = cholupdate(R, [2/sqrt(det); zeros(n-1,1)], '-');
                if err > 0
                    in = false;
                else
                    L = R.';
                    Li = @(M)(L\M);
                end
            end
        end
    end
    
return
