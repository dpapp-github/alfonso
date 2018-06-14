% This script demonstrates how to use the provided methods to solve the
% polynomial optimization problems described in:
%
% D. Papp and S. Yildiz. Sum-of-squares optimization without semidefinite 
% programming. Available at https://arxiv.org/abs/1712.01792.
% -------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz    <syildiz@email.unc.edu>  
%
% Date: 06/14/2018
%
% This code has been developed and tested with Matlab R2016b.
% -------------------------------------------------------------------------

tol = 1e-06;

intParams = PaduaSquare(10); % in the case of bivariate polynomials
% intParams = FeketeCube(4, 2); % in the case of n-variate polynomials

results = polyOpt(intParams, 'robinson', tol);
% results = polyOpt(intParams, 'caprasse', tol);
