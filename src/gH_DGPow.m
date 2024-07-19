%GH_DGPOW implements a membership and barrier function oracle for the dual
% of the generalized power cone
%   P^*_{a} = {(x,z)\in R_+^n x R | prod((x/a).^a) >= abs(z)},
% also known as the SONC (sum-of-nonnegative-circuit-polynomials cone)
%
% It is parametrized by the vector a which must satisfy sum(a) = 1 and
% a(:) > 0.
%
% --------------------------------------------------------------------------
% USAGE of "gH_DGPow"
% [in, g, Hi, Li] = gH_DGPow(xz, a)
% --------------------------------------------------------------------------
% INPUT
% xz:           primal iterate xz = [x; z]
% a:            parameter vector satisfying sum(a) = 1 and a(:) > 0
%
% OUTPUT
% in:	0 if x is not in the interior of the cone. 1 if x is in the
%       interior of the cone.
% g:	gradient of the barrier function at x
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
% Copyright (C) 2020 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz
% --------------------------------------------------------------------------


function [in, g, Hi, Li] = gH_DGPow(xz, a)

    xz = xz ./ [a;1];    % only for the dual cone
    x  = xz(1:end-1);
    z  = xz(end);

    %alogx = a'*log(x);
    %in    = all(x>0) & alogx > log(abs(z));
    xpa    = prod(x.^a);
    in     = all(x>0) & xpa > abs(z);
    
    g = NaN; Hi = NaN; Li = NaN; % will be overwritten if in==1

    if in
        %xpa    = exp(alogx);
        x2amz2 = xpa^2-z^2; %(xpa+z)*(xpa-z);

        % gradient; simple formula
        g = [-(1-a)./x - 2*(xpa^2.*a./x)/x2amz2; 2*z/x2amz2];
        g = g./[a;1];          % only for the dual cone

        % Hessian
        if nargout > 2
            % actual H
            H = 4*x2amz2^(-2)*xpa^2*([z*a;-1]*[z*a;-1]')./([x;1]*[x;1]') + ...
                diag( [(z^2*(-1+a)+xpa^2*(1+a))./(x2amz2*x.^2) ; -2/x2amz2] );
            H = H./([a;1]*[a;1]'); % only for the dual cone

            [L,err] = chol(H,'lower');
            if err > 0
                in = false;
            else
                Hi = @(v)( L'\(L\v) );
                Li = @(M)( L\M );
            end

        end % if nargout > 2
    end % if in

return

