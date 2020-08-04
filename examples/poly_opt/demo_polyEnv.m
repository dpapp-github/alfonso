% This script demonstrates how to use the provided methods to solve the
% polynomial envelope problem described in:
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

seed = 2017;

tol = 1e-06;

intParams = ChebInterval(100); % in the case of univariate polynomials
intParams = PaduaSquare(10); % in the case of bivariate polynomials
% intParams = FeketeCube(3, 6); % in the case of n-variate polynomials

results = polyEnv(intParams, 2, 5, tol, seed);
