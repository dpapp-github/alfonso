% This script demonstrates how to use the provided methods to solve a 
% random linear programming problem.
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

results = random_lp(500, 1000, tol, seed);
