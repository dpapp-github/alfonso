%GH_EXP implements a membership and barrier function oracle for the
% exponential cone
%   K_{exp} :=  cl { x \in R_+^2 x R | x1 > x2*exp(x3/x2) }
%
% It requires no parameters and it is always three-dimensional.
% --------------------------------------------------------------------------
% USAGE of "gH_Exp"
% [in, g, Hi, Li] = gH_Exp(x)
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
% None.
% -------------------------------------------------------------------------
% Copyright (C) 2020 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz
% --------------------------------------------------------------------------

function [in, g, Hi, Li] = gH_Exp(x, ~)

    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    
    logx1px2 = log(x1/x2);
    den      = -x3 + x2*logx1px2;

    g = NaN; Hi = NaN; Li = NaN; % will be overwritten if in==1
    
    in = x1 > 0 && x2 > 0 && den > 0;
    
    if in
        g = [-(1 + x2/den)/x1; (1-x3/x2)/den - 2/x2; 1/den];
 
        if nargout > 2
            % not really H yet
            H = [ (x2^2 + x2*den + den^2)/(x1^2), (-x2+x3)/x1, -x2/x1;
                  (-x2+x3)/x1, ((x2-x3)^2 + 2*den^2- (x2-2*x3)*den)/(x2^2) , 1-logx1px2;
                  -x2/x1, 1-logx1px2, 1];
        end
        
        if nargout >= 3
            [L, err] = chol(sparse(H), 'lower');
            
            if err > 0
                in = false;
                return;
            else
                %H = den^(-2)*sparse(H);
                %L = den^(-1)*L;
                Hi = @(v)( den^2 * (H\v) );
                Li = @(M)( den * (L\M) );
            end      
        end
    end
    
return
