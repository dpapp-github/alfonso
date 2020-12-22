% This code generates parameters for the interpolant basis representation
% of univariate sum-of-squares polynomials.
% -------------------------------------------------------------------------
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
% Version: 06/14/2018
%
% This code has been developed and tested with Matlab R2016b.
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% chebpts, chebpolyval from Chebfun. Chebfun is an open-source package for  
% numerical computing with functions: http://www.chebfun.org/.
% Our code has been developed and tested with Chebfun version 5.5.0.  
% The latest version of the Chebfun package can be downloaded from
% http://www.chebfun.org/download/.
% -------------------------------------------------------------------------


function intParams = ChebInterval(d)
% This method generates parameters for the interpolant basis representation
% of univariate sum-of-squares polynomials.
% --------------------------------------------------------------------------
% USAGE of "ChebInterval"
% intParams = ChebInterval(d)
% --------------------------------------------------------------------------
% INPUT
% d:	degree of polynomials to be squared
%
% OUTPUT
% intParams:        interpolation parameters
% - intParams.n:    number of arguments to the polynomials. intParams.n = 1
%                   in the univariate case.
% - intParams.d:    degree of polynomials to be squared
% - intParams.L:    dimension of the space of (intParams.n)-variate
%                   degree-d polynomials. intParams.L = d+1 in the
%                   univariate case.
% - intParams.U:    dimension of the space of (intParams.n)-variate
%                   degree-(2*d) polynomials. intParams.U = 2*d+1 in the
%                   univariate case.
% - intParams.pts:  Chebyshev points of the second kind for degree-(2*d)
%                   polynomial interpolation. (intParams.U x 1) array.
% - intParams.w:    (scaled) weights for Clenshaw-Curtis quadrature
% - intParams.P0:   evaluations of Chebyshev polynomials of the first kind
%                   up to degree d at the points intParams.pts. 
%                   (intParams.U x intParams.L) array.
% - intParams.P:    evaluations of a basis for the space of
%                   (intParams.n)-variate degree-d polynomials at the points
%                   intParams.pts. (intParams.U x intParams.L) array with
%                   orthonormal columns.
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% chebpts, chebpolyval from Chebfun. Chebfun is an open-source package for  
% numerical computing with functions: http://www.chebfun.org/.
% Our code has been developed and tested with Chebfun version 5.5.0.  
% The latest version of the Chebfun package can be downloaded from
% http://www.chebfun.org/download/.
% -------------------------------------------------------------------------
    
    n = 1; % univariate polynomials
    
    intParams.n = n;
    intParams.d = d;    
    intParams.L = nchoosek(n+d, n);
    intParams.U = nchoosek(n+2*d, n);
    [intParams.pts, intParams.w] = chebpts(intParams.U);
    intParams.w = intParams.w(:);
    intParams.P = zeros(intParams.U, intParams.L);
    
    col     = 0;
    lrEye   = fliplr(eye(d+1));
    for t = 0:d % polynomials with degree up to d
        col                 = col+1;
        intParams.P(:,col)  = chebpolyval(lrEye(:,t+1), intParams.pts);
    end
    
    intParams.P0 = intParams.P;
    [intParams.P, ~] = qr(intParams.P, 0);
    
return
    