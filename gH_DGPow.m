% gH_DGPow implements the usual logarithmic barrier function
% for the nonnegative circuit polynomial cone, aka. dual of the generalized power cone

%GH_DGPOW implements a membership and barrier function oracle for the dual
% of the generalized power cone
%   P^*_{a} = {(x,z)\in R_+^n x R | prod((x/a).^a) >= abs(z)},
% also known as the SONC (sum-of-nonnegative-circuit-polynomials cone)
%
% It is parametrized by the vector a which must satisfy sum(a) = 1 and
% a(:) > 0.
% --------------------------------------------------------------------------
% USAGE of "gH_DGPow"
% [in, g, H, L] = gH_DGPow(xz, a)
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

function [in, g, H, L] = gH_DGPow(xz, a)

    xz = xz ./ [a;1];    % only for the dual cone
    x  = xz(1:end-1);
    z  = xz(end);

    %alogx = a'*log(x);
    %in    = all(x>0) & alogx > log(abs(z));
    xpa    = prod(x.^a);
    in     = all(x>0) & xpa > abs(z);
    
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
            
            if nargout > 3
                [L,err] = chol(H,'lower');
                if err > 0
                    in = false; g = NaN; H = NaN; L = NaN;
                    return;
                end
            end
        end
    else
        g = NaN;
        H = NaN;
        L = NaN;
    end

return

