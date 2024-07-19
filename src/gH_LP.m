%GH_LP implements a membership and barrier function oracle for the
% nonnegative orthant, using the standard logarithmic barrier.
% --------------------------------------------------------------------------
% USAGE of "gH_LP"
% [in, g, Hi, Li] = gH_LP(x)
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
%
% --------------------------------------------------------------------------

function [in, g, Hi, Li] = gH_LP(x, ~)

    in = min(x)>0;

    if in
        % The gradient in closed form.
        g = -1./x;

        % Hessian-related functions in closed form.
        if nargout >= 3
            Hi = @(v)( (x.^2).*v );
            if nargout == 4
                Li = @(v)( v .* x );
            end          
        end
    else
        g = NaN; Li = NaN; Hi = NaN;
    end

return

