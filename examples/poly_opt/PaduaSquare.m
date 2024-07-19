% This code generates parameters for the interpolant basis representation
% of bivariate sum-of-squares polynomials.
% -------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz
%
% Date: 06/14/2018
%
% This code has been developed and tested with Matlab R2023b.
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% pdpts, pdwtsMM from Padua2DM. Padua2DM is a Matlab package from M. Caliari,
% S. De Marchi, A. Sommariva, and M. Vianello for interpolation and
% cubature at the Padua points. It can be downloaded from
% http://profs.sci.univr.it/~caliari/software.htm.
%
% chebpolyval from Chebfun. Chebfun is an open-source package for  
% numerical computing with functions: http://www.chebfun.org/.
% Our code has been developed and tested with Chebfun version 5.5.0.  
% The latest version of the Chebfun package can be downloaded from
% http://www.chebfun.org/download/.
%
% partitions. The partitions function computes all partitions of an integer.
% We use the implementation of John D'Errico from a MATLAB Central File 
% Exchange post which is available at
% https://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer.
% -------------------------------------------------------------------------


function intParams = PaduaSquare(d)
% This method generates parameters for the interpolant basis representation
% of bivariate sum-of-squares polynomials.
% --------------------------------------------------------------------------
% USAGE of "PaduaSquare"
% intParams = PaduaSquare(d)
% --------------------------------------------------------------------------
% INPUT
% d:    degree of polynomials to be squared
%
% OUTPUT
% intParams:        interpolation parameters
% - intParams.n:    number of arguments to the polynomials. intParams.n = 2
%                   in the bivariate case.
% - intParams.d:    degree of polynomials to be squared
% - intParams.L:    dimension of the space of (intParams.n)-variate
%                   degree-d polynomials. intParams.L = nchoosek(d+2,2) in
%                   the bivariate case.
% - intParams.U:    dimension of the space of (intParams.n)-variate
%                   degree-(2*d) polynomials. intParams.U = nchoosek(2*d+2,2)
%                   in the bivariate case.
% - intParams.pts:  Padua points for degree-(2*d) polynomial interpolation. 
%                   (intParams.U x 1) array.
% - intParams.w:    (scaled) weights for Clenshaw-Curtis quadrature
% - intParams.P0:   evaluations of bivariate product Chebyshev polynomials
%                   of the first kind up to degree d at the points
%                   intParams.pts. (intParams.U x intParams.L) array.
% - intParams.P:    evaluations of a basis for the space of
%                   (intParams.n)-variate degree-d polynomials at the points
%                   intParams.pts. (intParams.U x intParams.L) array with
%                   orthonormal columns.
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% pdpts, pdwtsMM from Padua2DM. Padua2DM is a Matlab package from M. Caliari,
% S. De Marchi, A. Sommariva, and M. Vianello for 
% interpolation and cubature at the Padua points. It can be downloaded from
% http://profs.sci.univr.it/~caliari/software.htm.
%
% chebpolyval from Chebfun. Chebfun is an open-source package for  
% numerical computing with functions: http://www.chebfun.org/.
% Our code has been developed and tested with Chebfun version 5.5.0.  
% The latest version of the Chebfun package can be downloaded from
% http://www.chebfun.org/download/.
%
% partitions. The partitions function computes all partitions of an integer.
% We use the implementation of John D'Errico from a MATLAB Central File 
% Exchange post which is available at
% https://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer.
% -------------------------------------------------------------------------

    n = 2;  % bivariate polynomials
    
    intParams.n = n;
    intParams.d = d;    
    intParams.L = nchoosek(n+d,n);
    intParams.U = nchoosek(n+2*d,n);
    intParams.pts = pdpts(2*d);
    intParams.w = pdwtsMM(2*d);
    intParams.P = ones(intParams.U,intParams.L);
        
    col = 0;
    lrEye = fliplr(eye(d+1));
    for t = 0:d     % polynomials with total degree up to d
        allDegs = partitions(t, ones(1,n));
        [nrDegs,~] = size(allDegs);
        for i = 1:nrDegs
            col = col+1;
            for j = 1:n
                dj = allDegs(i,j);
                intParams.P(:,col) = ...
                    intParams.P(:,col).*chebpolyval(lrEye(:,dj+1),intParams.pts(:,j));
            end
        end
    end
    
    intParams.P0 = intParams.P;
    [intParams.P,~] = qr(intParams.P,0);
    
return
