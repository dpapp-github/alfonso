% This code formulates and solves a random linear programming problem.
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

function results = random_lp(m, n, tol, seed)
% This is the main method for solving a random linear programming problem.
% --------------------------------------------------------------------------
% USAGE of "random_lp"
% results = random_lp(m, n, tol, seed)
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


    rng(seed, 'twister');

    A  = randi([-9,9],m,n);
    b  = A*ones(n, 1);
    c  = randi(9,n,1);

    probData = struct('c', c, 'A', A, 'b', b);

    tic
    
    % INITIAL PRIMAL ITERATE
    x0 = ones(n,1);
    [~, g0] = gH_LP(x0, []);
    % scaling factor for the primal problem
    rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
    % scaling factor for the dual problem
    rD = max((1+abs(g0))./(1+abs(probData.c)));
    % initial primal iterate
    x0 = repmat(sqrt(rP*rD),n,1);
    
    % CUSTOM ALGORITHMIC OPTIONS
    opts.optimTol = tol;
        
    % CALL TO alfonso 
    results = alfonso(probData, x0, @gH_LP, [], opts);
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
