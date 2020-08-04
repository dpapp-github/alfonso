% This code formulates and solves a random linear programming problem
% utilizing the simple interface.
% -------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz    <syildiz@email.unc.edu>  
%
% Date: 01/15/2019
%
% This code has been developed and tested with Matlab R2016b.
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% None.
% -------------------------------------------------------------------------

function results = random_lp_simple(m, n, tol, seed)
% This is the main method for solving a random linear programming problem.
% --------------------------------------------------------------------------
% USAGE of "random_lp_simple"
% results = random_lp_simple(m, n, tol, seed)
% --------------------------------------------------------------------------
% INPUT
% m:	number of equality constraints
% n:	number of variables
% tol:	tolerance parameter for optimization
% seed:	seed for the random number generator
%
% OUTPUT
% results:                  final solution and iteration statistics
% - results.nIterations:	total number of iterations
% - results.alphaPred:      predictor step size at each iteration
% - results.betaPred:       neighborhood parameter at the end of the
%                           predictor phase at each iteration
% - results.etaCorr:        neighborhood parameter at the end of the
%                           corrector phase at each iteration
% - results.mu:             complementarity gap at each iteration
% - results.x:              final value of the primal variables
% - results.s:              final value of the dual slack variables
% - results.y:              final value of the dual free variables
% - results.tau:            final value of the tau-variable
% - results.kappa:          final value of the kappa-variable
% - results.pObj:           final primal objective value
% - results.dObj:           final dual objective value
% - results.dGap:           final duality gap
% - results.cGap:           final complementarity gap
% - results.rel_dGap:       final relative duality gap   
% - results.rel_cGap:       final relative complementarity gap
% - results.pRes:           final primal residuals
% - results.dRes:           final dual residuals
% - results.pIn:            final primal infeasibility
% - results.dIn:            final dual infeasibility
% - results.rel_pIn:        final relative primal infeasibility
% - results.rel_dIn:        final relative dual infeasibility
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    if isOctave()
        rand("state", seed);
    else
        rng(seed, 'twister');
    end

    A  = randi([-9,9],m,n);
    b  = A*ones(n,1);
    c  = randi(9,n,1);

    tic
    
    % INITIAL PRIMAL ITERATE
    x0 = ones(n,1);
    [~, g0] = gH_lp(x0, []);
    % scaling factor for the primal problem
    rP = max((1+abs(b))./(1+abs(A*x0)));
    % scaling factor for the dual problem
    rD = max((1+abs(g0))./(1+abs(c)));
    % initial primal iterate
    x0 = repmat(sqrt(rP*rD),n,1);
    
    % CUSTOM ALGORITHMIC OPTIONS
    opts = struct('optimTol', tol, 'preprocess', false);
        
    % CALL TO alfonso 
    K{1} = struct('type', 'lp', 'dim', n);
    results = alfonso_simple(c, A, b, K, x0, opts);
    fprintf('alfonso is done. ');
    
    toc
    
    % prints error statistics
    fprintf('\n');
    fprintf('FINAL:\n');
    fprintf('Relative primal infeasibility: %d\n', results.rel_pIn);
    fprintf('Relative dual infeasibility: %d\n', results.rel_dIn);
    fprintf('Relative duality gap: %d\n', results.rel_dGap);
    fprintf('Relative complementarity gap: %d\n\n', results.rel_cGap);

end

function [in, g, H, L] = gH_lp(x, ~)
% This method computes the gradient and Hessian of the barrier function for
% the linear program. It requires no parameters.
% --------------------------------------------------------------------------
% USAGE of "gH_lp"
% [in, g, H, L] = gH_lp(x)
% --------------------------------------------------------------------------
% INPUT
% x:            primal iterate
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

    n  = length(x);
    in = min(x)>0;
    
    if in
        g = -1./x;
        H = sparse(1:n,1:n,x.^(-2),n,n,n);
        L = sparse(1:n,1:n,-g,n,n,n);
    else
        g = NaN; H = NaN; L = NaN;
    end

end
