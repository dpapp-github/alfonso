%GH_LP implements a membership and barrier function oracle for the
% nonnegative orthant, using the standard logarithmic barrier. It requires
% no parameters.
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
% USAGE of "gH_LP"
% [in, g, H, L] = gH_LP(x)
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
