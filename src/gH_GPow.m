%GH_GPOW implements a membership and barrier function oracle for the
% generalized power cone
%   P_{a} = {(x,z)\in R_+^n x R | prod(x.^a) >= abs(z)}
%
% It is parametrized by the vector a which must satisfy \sum(a) = 1 and
% a > 0.
% --------------------------------------------------------------------------
% USAGE of "gH_GPow"
% [in, g, Hi, Li] = gH_GPow(xz, a)
% --------------------------------------------------------------------------
% INPUT
% xz:       primal iterate xz = [x; z]
% a:        parameter vector satisfying sum(a) = 1 and a(:) > 0
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

function [in,g,Hi,Li] = gH_GPow(xz, a)

    x  = xz(1:end-1);
    z  = xz(end);

    g = NaN; Hi = NaN; Li = NaN;  % will be overwritten if in==1

    xpa = prod(x.^a);
    in  = all(x>0) & xpa > abs(z);

    if in
        x2amz2 = xpa^2-z^2; %(xpa+z)*(xpa-z);

        % gradient; simple formula
        g = [-(1-a)./x - 2*(xpa^2.*a./x)/x2amz2; 2*z/x2amz2];

        % Hessian
        if nargout >= 3
            H = 4*x2amz2^(-2)*xpa^2*([z*a;-1]*[z*a;-1]')./([x;1]*[x;1]') + ...
                diag( [(z^2*(-1+a)+xpa^2*(1+a))./(x2amz2*x.^2) ; -2/x2amz2] );
            
            % direct solve with Cholesky
            [L,err] = chol(H,'lower');
            if err
                in = 0;
            else
                Hi = @(v)(L'\(L\v));
                if nargout == 4
                    Li = @(M) L\M;
                end % if nargout == 4
            end
        end
    end

return
