% This script demonstrates how to use the provided methods to solve the
% polynomial optimization problems described in:
%
% D. Papp and S. Yildiz. Sum-of-squares optimization without semidefinite 
% programming. SIAM Journal on Optimization 29(1), 2019, pp. 822-851.
% URL: https://doi.org/10.1137/17M1160124, preprint available at
% https://arxiv.org/abs/1712.01792. 
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

tol = 1e-06;

% intParams = PaduaSquare(10);    % in the case of bivariate polynomials
intParams = FeketeCube(4, 2); % in the case of n-variate polynomials; the
                                % first argument must be n

% results = polyOpt(intParams, 'robinson', tol);   % bivariate example
results = polyOpt(intParams, 'caprasse', tol);     % 4-variate example
