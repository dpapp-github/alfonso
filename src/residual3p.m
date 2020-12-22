% This code is an implementation of the triple-precision accumulated inner
% product and residual computation. It was adapted from the code of 
% Cleve Moler from a MATLAB Central Blog post which is available at
% https://blogs.mathworks.com/cleve/2015/03/02/triple-precision-accumlated-inner-product/
% -------------------------------------------------------------------------
% Version: 06/14/2018
%
% This code has been developed and tested with Matlab R2016b.
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% None.
% --------------------------------------------------------------------------


function res = residual3p(LHS,delta,RHS)
% This method computes the residual "LHS*delta - RHS" in triple precision.
% --------------------------------------------------------------------------
% USAGE of "residual3p"
% r = residual3p(LHS,delta,RHS)
% --------------------------------------------------------------------------
% INPUT
% LHS:     left-hand side matrix
% delta:   solution vector
% RHS:     right-hand side vector
%
% OUTPUT
% res:     residual vector
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    m = size(LHS,1);
    res = zeros(m,1);   

    for k = 1:m
        res(k,:) = dot3p(LHS(k,:),delta,-RHS(k,:));
    end

return

function total = dot3p(x,y,total)
% This method computes the extended inner product "x*y + total" in 
% triple precision.
% --------------------------------------------------------------------------
% USAGE of "dot3p"
% total = dot3p(x,y,total)
% --------------------------------------------------------------------------
% INPUT
% x:        row vector
% y:        column vector
% total:    scalar
%
% OUTPUT
% total:    result of the inner product
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    % Matlab does not allow sparse matrices at single precision
    if issparse(x); x = full(x); end;
    
    shi = double(single(total));
    slo = total - shi;
    for k = 1:length(x)
        xhi = double(single(x(k)));
        xlo = x(k) - xhi;
        yhi = double(single(y(k)));
        ylo = y(k) - yhi;
        tmp = xhi*yhi;
        zhi = double(single(tmp));
        zlo = tmp - zhi + xhi*ylo + xlo*yhi + xlo*ylo;

        tmp = shi + zhi;
        del = tmp - shi - zhi;
        shi = double(single(tmp));
        slo = tmp - shi + slo + zlo - del;
    end
    total = shi +  slo;
   
return
