%GH_GPOW implements a membership and barrier function oracle for the
% generalized power cone
%   P_{a} = {(x,z)\in R_+^n x R | prod(x.^a) >= abs(z)}
%
% It is parametrized by the vector a which must satisfy \sum(a) = 1 and
% a > 0.
% --------------------------------------------------------------------------
% USAGE of "gH_GPow"
% [in, g, H, L] = gH_GPow(xz, a)
% --------------------------------------------------------------------------
% INPUT
% xz:           primal iterate xz = [x; z]
% a:            parameter vector satisfying sum(a) = 1 and a(:) > 0
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

function [in,g,H,L] = gH_GPow(xz, a)

    x  = xz(1:end-1);
    z  = xz(end);

    xpa    = prod(x.^a);
    in     = all(x>0) & xpa > abs(z);

    if in
        x2amz2 = xpa^2-z^2; %(xpa+z)*(xpa-z);

        % gradient; simple formula
        g = [-(1-a)./x - 2*(xpa^2.*a./x)/x2amz2; 2*z/x2amz2];

        % Hessian
        if nargout > 2
            H = 4*x2amz2^(-2)*xpa^2*([z*a;-1]*[z*a;-1]')./([x;1]*[x;1]') + ...
                diag( [(z^2*(-1+a)+xpa^2*(1+a))./(x2amz2*x.^2) ; -2/x2amz2] );
        end

        if nargout > 3
            [L,err] = chol(H,'lower');
            if err > 0
                in = false; g = NaN; H = NaN; L = NaN;
                return;
            end
        end
    else
        g = NaN;
        H = NaN;
        L = NaN;
    end

return
