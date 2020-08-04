% This code is an implementation of the algorithm for non-symmetric conic 
% optimization, which was originally presented in:
%
% A. Skajaa and Y. Ye, A homogeneous interior-point algorithm for nonsymmetric 
% convex conic optimization, Mathematical Programming Ser. A, 150 (2015), 
% pp. 391-422. Available at https://doi.org/10.1007/s10107-014-0773-1.
%
% The implementation is based on the corrected analysis of the algorithm
% presented in:
%
% D. Papp and S. Yildiz. On "A homogeneous interior-point algorithm for
% nonsymmetric convex conic optimization". Available at 
% https://arxiv.org/abs/1712.00492.
% -------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz    <syildiz@email.unc.edu>  
%
% Version: 2020/07/20
% 
% This code has been developed and tested with Matlab R2016b and
% Octave v.4.4.1.
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% None.
% --------------------------------------------------------------------------


function results = alfonso(probData, x0, gH, gH_Params, opts)
% ALgorithm FOr Non-Symmetric Optimization
% This is the main method for the algorithm.
% --------------------------------------------------------------------------
% USAGE of "alfonso"
% results = alfonso(probData, x0, gH, gH_Params, opts)
% --------------------------------------------------------------------------
% INPUT
% probData:                     data for the conic optimization problem
% - probData.A:                 constraint matrix
% - probData.b:                 right-hand side vector
% - probData.c:                 cost vector
% x0:                           initial primal iterate
% gH:                           method for computing the gradient and 
%                               the factored Hessian of the barrier function
% gH_Params:                    parameters required by the method gH; may
%                               be replaced by [] if no parameters needed
% opts:                         algorithmic options
% - opts.maxIter:               maximum number of interior-point iterations.
%                               default value: 10000
% - opts.predLineSearch:        0 if a fixed step size is to be used in the 
%                               predictor step. 1 if the step size is to be 
%                               determined via line search in the predictor
%                               step. default value: 1.
% - opts.maxCorrSteps:          maximum number of corrector steps. 
%                               possible values: 1, 2, or 4. default value: 4.
% - opts.corrCheck:             0 if maxCorrSteps corrector steps are to be 
%                               performed at each corrector phase. 1 if the
%                               corrector phase can be terminated before 
%                               maxCorrSteps corrector steps if the iterate 
%                               is in the eta-neighborhood. default value: 1.
% - opts.optimTol:              optimization tolerance parameter. default
%                               value: 1e-06. minimum value: eps.
% - opts.maxCorrLSIters:        maximum number of line search iterations in
%                               each corrector step. default value: 8.
% - opts.maxSmallPredSteps:     maximum number of predictor step size 
%                               reductions allowed with respect to the safe
%                               fixed step size. default value: 8.
% - opts.maxItRefineSteps:      maximum number of iterative refinement steps 
%                               in linear system solves. default value: 0.
% - opts.verbose:               0 if output is to be suppressed. 1 if
%                               progress is to be printed at each iteration.
%                               default value: 1.
%
% OUTPUT
% results:                  final solution and iteration statistics
% - results.status:         solver status: 1 = success, 0 = infeasible problem, everything else is trouble
% - results.statusString:   solver status string
% - results.nIterations:	total number of iterations
% - results.x:              final value of the primal variables
% - results.s:              final value of the dual slack variables
% - results.y:              final value of the dual free variables
% - results.tau:            final value of the tau-variable
% - results.kappa:          final value of the kappa-variable
% - results.pObj:           final primal objective value
% - results.dObj:           final dual objective value
% - results.options:        options structure used in the computation
% - results.alphaPred:      predictor step size at each iteration
% - results.betaPred:       neighborhood parameter at the end of the
%                           predictor phase at each iteration
% - results.etaCorr:        neighborhood parameter at the end of the
%                           corrector phase at each iteration
% - results.mu:             complementarity gap at each iteration
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
% --------------------------------------------------------------------------

    if nargin == 4
        opts = struct();
    elseif nargin < 4
        error('alfonso needs more input arguments.')
    elseif nargin > 5
        error('alfonso needs fewer input arguments.')
    end
    
    % sets algorithmic options
    opts = setOpts(opts);
    
    % say hello, alfonso
    say_hello(opts);
    
    % checks the problem data for consistency
    inputCheck(probData);

    % check the initial point and compute the barrier parameter of the barrier function
    [in, g] = gH(x0, gH_Params);
    if ~in
        error('Specified initial point is not in the cone.');
    end
    bnu = (-g'*x0) + 1;  % nu-bar = nu+1, where nu = g(x0)'*x0 is the barrier parameter
    if isfield(gH_Params,'bnu')
        if abs(gH_Params.bnu-bnu) > 1e-12
            warning('Specified and computed gH_Params.bnu arguments do not agree. Specified: %d; computed: %d', gH_Params.bnu, bnu);
        end
    else
        gH_Params.bnu = bnu;
        if opts.verbose
            disp(['barrier parameter set to nu = ', num2str(bnu-1)]);
        end
    end
    
    
    [m, n] = size(probData.A);
    A = probData.A;
    b = probData.b;
    c = probData.c;
    
    if opts.maxItRefineSteps > 0 
        probData.LHS = ...
        [ sparse(m,m)   A               -b            sparse(m,n)    sparse(m,1) ;
         -A'            sparse(n,n)      c           -speye(n)       sparse(n,1) ;
          b'           -c'               0            sparse(1,n)   -1          ;
          sparse(n,m)   speye(n)         sparse(n,1)  speye(n)       sparse(n,1) ;
          sparse(1,m)   sparse(1,n)      1            sparse(1,n)    1          ];
    end
    
    % sets the solution method for the Newton system
    myLinSolve = @linSolve5;
   
    % sets algorithmic parameters
    algParams = setAlgParams(gH_Params, opts);
    
    results.status = 1;
    results.statusString = '';
    
    % creates arrays for iteration statistics
    results.alphaPred   = zeros(algParams.maxIter, 1);
    results.betaPred    = zeros(algParams.maxIter, 1);
    results.etaCorr     = zeros(algParams.maxIter, 1);
    results.mu          = zeros(algParams.maxIter, 1);

    % sets constants for termination criteria
    termConsts.pRes = max([1, norm([A,b],Inf)]);
    termConsts.dRes = max([1, norm([A',speye(n),-c],Inf)]);
    termConsts.comp = max([1, norm([-c',b',1],Inf)]);

    % creates the central primal-dual iterate corresponding to x0
    soln = initSoln(x0, probData, gH, gH_Params);
    
    if ~opts.debug
        if isOctave()
            warning('off','Octave:nearly-singular-matrix');
        else
            warning('off','MATLAB:nearlySingularMatrix');
        end
    end
    termFlag  = 0;
    numIters  = 0;
    stopwatch = tic;
    elapsed   = 0;
  
    for iter = 1:algParams.maxIter+1

        % checks progress towards termination criteria
        [status, metrics] = term(soln, probData, algParams, termConsts);
            
        if termFlag || iter == algParams.maxIter+1 || (status ~= -99 && status ~= -6) 
            if iter == algParams.maxIter+1
                status = 0; % 'Number of iterations exceeded opts.maxIter.';
            end
            
            numIters = iter-1;
            break;
        end
        
        % prints progress metrics
        if mod(iter,1)==0 && opts.verbose
            fprintf('%3d: pObj=%.6e pIn=%#.2e dIn=%#.2e gap=%#.2e tau=%#.2e kap=%#.2e mu=%.2e t=%#.2f s\n',...
            iter, metrics.O, metrics.P, metrics.D, metrics.A, soln.tau, soln.kappa, soln.mu, elapsed);
        end

        % PREDICTOR PHASE
        [soln, alphaPred, betaPred, algParams, predStatus] =...
            pred(soln, probData, gH, gH_Params, myLinSolve, algParams, opts);
        results.alphaPred(iter) = alphaPred;
        results.betaPred(iter)  = betaPred;
        % raises a termination flag if predictor phase was not successful
        if predStatus == 0
            if iter > 1
                results.betaPred(iter) = results.etaCorr(iter-1);
            end
            if opts.verbose
                fprintf('Predictor could not improve the solution.\n');
            end
            termFlag = 1;
        end

        % CORRECTOR PHASE
        results.etaCorr(iter) = results.betaPred(iter);
        % skips corrector phase if corrCheck == 1 and current iterate is 
        % in the eta-neighborhood 
        if (~opts.corrCheck || results.etaCorr(iter) > algParams.eta) && ~termFlag
            for corrIter = 1:algParams.maxCorrSteps
                [soln, corrStatus] = corr(soln, probData, gH, gH_Params, myLinSolve, algParams, opts);
                % exits corrector phase and raises a termination flag if 
                % last corrector step was not successful
                if corrStatus == 0
                    if opts.verbose
                        fprintf('Corrector could not improve the solution.\n');
                    end
                    termFlag = 1; break;
                end
                % exits corrector phase if corrCheck == 1 and current
                % iterate is in the eta-neighborhood
                if (opts.corrCheck && corrIter < algParams.maxCorrSteps) || corrIter == algParams.maxCorrSteps
                    results.etaCorr(iter) = sqrt(sum((soln.L\soln.psi(1:end-1)).^2) +...
                        (soln.tau*soln.psi(end))^2)/soln.mu;
                    if results.etaCorr(iter) <= algParams.eta; break; end
                end
            end
            % raises a termination flag if corrector phase was not successful
            if opts.debug && results.etaCorr(iter) > algParams.eta
                if opts.verbose
                    fprintf('Corrector phase finished outside the eta-neighborhood.\n');
                end
                %termFlag = 1;      
            end
        end
        results.mu(iter) = soln.mu;
        
        elapsed = toc(stopwatch);
    end
    
    if isOctave()
        warning('on','Octave:nearly-singular-matrix');
    else
        warning('on','MATLAB:nearlySingularMatrix');
    end
    % prepares final solution and iteration statistics
    results = prepResults(results, status, soln, probData, numIters, elapsed, opts);

    if opts.verbose
        disp(['Done in ', int2str(numIters), ' iterations.']);
        disp(['Status = ', results.statusString]);
    end
return

function [solnAlpha, alpha, betaAlpha, algParams, predStatus] = pred(soln, probData, gH, gH_Params, myLinSolve, algParams, opts)
% This method performs a predictor step in the algorithm.
% --------------------------------------------------------------------------
% USAGE of "pred"
% [solnAlpha, alpha, betaAlpha, algParams, predStatus] = pred(soln, probData, gH, gH_Params, myLinSolve, algParams, opts)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% gH:           method for computing the gradient and Hessian of the barrier function
% gH_Params:	parameters associated with the method gH
% myLinSolve:   solution method for the Newton system
% algParams:    algorithmic parameters
% opts:         algorithmic options
%
% OUTPUT
% solnAlpha:    new iterate
% alpha:        predictor step size
% betaAlpha:    neighborhood parameter at the end of the predictor phase
% algParams:    algorithmic parameters
% predStatus:   0 if predictor phase was not successful. 1 if predictor 
%               phase was successful.
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    x       = soln.x;
    tau     = soln.tau;
    s       = soln.s;
    kappa   = soln.kappa;
    y       = soln.y;
    
    A = probData.A;
    b = probData.b;
    c = probData.c;
    [m, n] = size(A);
    
    RHS = zeros(m+2*n+2,1);
    RHS(1:m+n+1)        = -[A*x - b*tau; -A'*y + c*tau - s; b'*y - c'*x - kappa];
    RHS(m+n+2:end)      = -[s; kappa];
    
    % computes Newton direction
    dsoln = linSolveMain(soln, probData, RHS, myLinSolve, algParams, opts);
    
    % line search is done even when predLineSearch == 0 to make sure next
    % iterate stays within the beta-neighborhood.
    betaAlpha       = Inf; 
    solnAlphaInNhd  = 0;
    alphaPrevOK     = -1;
    predStatus      = 1;
    
    nSteps = 0;
    alpha = algParams.alphaPred;
    
    while true
        nSteps = nSteps + 1;
        
        solnAlpha.y         = y     + alpha*dsoln.y;
        solnAlpha.x         = x     + alpha*dsoln.x;
        solnAlpha.tau       = tau   + alpha*dsoln.tau;
        solnAlpha.s         = s     + alpha*dsoln.s;
        solnAlpha.kappa     = kappa + alpha*dsoln.kappa;            

        [solnAlpha.in, solnAlpha.g, solnAlpha.H, solnAlpha.L] = gH(solnAlpha.x, gH_Params);
        
        % primal iterate is inside the cone
        if solnAlpha.in
            solnAlpha.mu	= (solnAlpha.x'*solnAlpha.s +...
                solnAlpha.tau*solnAlpha.kappa)/gH_Params.bnu;
            solnAlpha.psi	= [solnAlpha.s; solnAlpha.kappa] +...
                solnAlpha.mu*[solnAlpha.g; -1/solnAlpha.tau];                
            betaAlpha       = sqrt(sum((solnAlpha.L\solnAlpha.psi(1:end-1)).^2) + ...
                (solnAlpha.tau*solnAlpha.psi(end))^2)/solnAlpha.mu;
            solnAlphaInNhd  = (betaAlpha < algParams.beta);
        end
        
        % iterate is inside the beta-neighborhood
        if solnAlpha.in && solnAlphaInNhd 
            % either the previous iterate was outside the beta-neighborhood
            % or increasing alpha again will make it > 1 
            if alphaPrevOK == 0 || alpha > algParams.predLSMulti
                if opts.predLineSearch == 1; algParams.alphaPred = alpha; end
                break;
            end
            alphaPrevOK     = 1;
            alphaPrev       = alpha;
            betaAlphaPrev   = betaAlpha;
            solnAlphaPrev   = solnAlpha;
            alpha           = alpha/algParams.predLSMulti;
        else    % iterate is outside the beta-neighborhood
            % previous iterate was in the beta-neighborhood
            if alphaPrevOK == 1
                alpha       = alphaPrev;
                betaAlpha   = betaAlphaPrev;
                solnAlpha   = solnAlphaPrev;
                if opts.predLineSearch == 1; algParams.alphaPred = alpha; end
                break;
            end
            % last two iterates were outside the beta-neighborhood and
            % alpha is very small
            if alpha < algParams.alphaPredThreshold
                predStatus  = 0; % predictor has failed
                alpha       = 0;
                betaAlpha   = Inf; % overwritten later in alfonso
                solnAlpha   = soln;
                if opts.predLineSearch == 1; algParams.alphaPred = alpha; end
                break;
            end
            % alphaPrev, betaAlphaPrev, solnAlphaPrev will not be used
            % given alphaPrevOK == 0
            alphaPrevOK     = 0;
            alphaPrev       = alpha;
            betaAlphaPrev   = betaAlpha;
            solnAlphaPrev   = solnAlpha;
            alpha = algParams.predLSMulti*alpha;
        end
    end
        
return

function [solnAlpha, corrStatus] = corr(soln, probData, gH, gH_Params, myLinSolve, algParams, opts)
% This method performs a single corrector step in the algorithm.
% --------------------------------------------------------------------------
% USAGE of "corr"
% [solnAlpha, corrStatus] = corr(soln, probData, gH, gH_Params, myLinSolve, algParams, opts)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% gH:           method for computing the gradient and Hessian 
%               of the barrier function
% gH_Params:    parameters associated with the method gH
% myLinSolve:   solution method for the Newton system
% algParams:    algorithmic parameters
% opts:         algorithmic options
%
% OUTPUT
% solnAlpha:    new iterate
% corrStatus:   0 if corrector phase was not successful. 1 if corrector 
%               phase was successful.
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    [m, n] = size(probData.A);

    RHS = zeros(m+2*n+2,1);
    RHS(m+n+2:end)  = -soln.psi;
    
    dsoln = linSolveMain(soln, probData, RHS, myLinSolve, algParams, opts);
    alpha = algParams.alphaCorr;
    corrStatus = 1;
    
    % does line search to make sure next primal iterate remains inside the cone
    for nSteps = 1:opts.maxCorrLSIters
        
        solnAlpha.y     = soln.y     + alpha*dsoln.y;
        solnAlpha.x     = soln.x     + alpha*dsoln.x;
        solnAlpha.tau   = soln.tau   + alpha*dsoln.tau;
        solnAlpha.s     = soln.s     + alpha*dsoln.s;
        solnAlpha.kappa = soln.kappa + alpha*dsoln.kappa;
        
        [solnAlpha.in, solnAlpha.g, solnAlpha.H, solnAlpha.L] = gH(solnAlpha.x, gH_Params);
        
        % terminates line search if primal iterate is inside the cone
        if solnAlpha.in
            solnAlpha.mu   = (solnAlpha.x'*solnAlpha.s + solnAlpha.tau*solnAlpha.kappa)/gH_Params.bnu;
            solnAlpha.psi  = [solnAlpha.s;solnAlpha.kappa] + solnAlpha.mu*[solnAlpha.g;-1/solnAlpha.tau];                
            break;
        end
        
        alpha = algParams.corrLSMulti*alpha;        
    end
        
    if solnAlpha.in == 0
        corrStatus  = 0; % corrector has failed
        solnAlpha   = soln;
    end
    
return

function dsoln = linSolveMain(soln, probData, RHS, myLinSolve, algParams, opts)
% This method sets up the Newton system and computes its solution.
% --------------------------------------------------------------------------
% USAGE of "linSolveMain"
% dsoln = linSolveMain(soln, probData, RHS, myLinSolve, algParams, opts)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% RHS:          right-hand side of the Newton system
% myLinSolve:   solution method for the Newton system
% algParams:    algorithmic parameters
% opts:         algorithmic options
%
% OUTPUT
% dsoln:	computed Newton direction
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------
    
    [m, n] = size(probData.A);
    delta  = myLinSolve(soln, probData, RHS);
    
    if opts.maxItRefineSteps > 0        
               
        % checks to see if we need to refine the solution
        if rcond(full(soln.H)) < eps            
            
            LHS                         = probData.LHS;    
            LHS(m+n+1+(1:n),m+(1:n))    = soln.mu*soln.H;
            %LHS(m+n+1+(1:n),m+(1:n))    = soln.mu*(soln.L*soln.L');
            LHS(end,m+n+1)              = soln.mu/soln.tau^2;
            
            exitFlag    = 0;
            res         = residual3p(LHS, delta, RHS);
            resNorm     = norm(res);  
            for iter = 1:opts.maxItRefineSteps
                if exitFlag; break; end
                d           = myLinSolve(soln, probData, res);
                deltaNew	= delta - d;
                resNew      = residual3p(LHS, deltaNew, RHS);
                resNewNorm	= norm(resNew);
                    
                % stops iterative refinement if there is not enough progress
                if resNewNorm > algParams.itRefineThreshold*resNorm
                    exitFlag = 1;
                end
                % updates solution if residual norm is smaller
                if resNewNorm < resNorm
                    delta       = deltaNew;
                    res         = resNew;
                    resNorm     = resNewNorm;
                end
            end
        end
        
    end
    
    dsoln.y     = delta(1:m);
    dsoln.x     = delta(m+(1:n));
    dsoln.tau   = delta(m+n+1);
    dsoln.s     = delta(m+n+1+(1:n));
    dsoln.kappa = delta(end);
    
return

function [delta, probData] = linSolve5(soln, probData, RHS)
% This method implements a block solver for the Newton system.
% --------------------------------------------------------------------------
% USAGE of "linSolve5"
% delta = linSolve5(soln, probData, RHS)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% RHS:          right-hand side of the Newton system
%
% OUTPUT
% delta:	computed Newton direction
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------
       
    A = probData.A;
    b = probData.b;
    c = probData.c;
    [m, n] = size(A);
    
    if issparse(A) && issparse(soln.L)
        L = sparse(soln.L);
        b = sparse(b);
        c = sparse(c);
        RHS = sparse(RHS);
    else
        L = full(soln.L);
    end
    
    mu     = soln.mu;
    tau    = soln.tau;
    
    ry     = RHS(1:m);
    rx     = RHS(m+(1:n));
    rtau   = RHS(m+n+1);
    rs     = RHS(m+n+1+(1:n));
    rkappa = RHS(end);

    LiAt    = L\A';
    Lic     = L\c;
    AHic    = LiAt'*Lic;
    
    Hirxrs      = L'\(L\(rx+rs));
    RHSdydtau   = [ry; rtau+rkappa] - [A; -c']*Hirxrs/soln.mu;
    
    x = RHSdydtau(1:end-1);
    y = RHSdydtau(end);
    
    dinv = (mu/tau^2 + Lic'*Lic/mu)^(-1);
    if isOctave() || nnz(LiAt)/numel(LiAt) > 0.1
        dy = (LiAt'*LiAt/mu + (b+AHic/mu)*dinv*(b-AHic/mu)') \ (x+(b+AHic/mu)*dinv*y);
    else
        dy = rk1upsolve(LiAt/sqrt(mu), -(b+AHic/mu)*dinv, b-AHic/mu, (x+(b+AHic/mu)*dinv*y));
    end
    dtau = dinv*((-b + AHic/mu)'*dy + y);
    
    dx = (Hirxrs + L'\([LiAt, -Lic]*[dy; dtau]))/soln.mu;
    
    delta               = zeros(m+2*n+2, 1);
    delta(1:m)          = dy;
    delta(m+n+1)        = dtau;
    delta(m+(1:n))      = dx;
    delta(m+n+1+(1:n))  = -rx - [A', -c]*[dy; dtau];
    delta(end)          = -rtau + b'*dy - c'*dx;
    
return

function x = rk1upsolve(Z, u, v, b)
% subroutine for linsolve5
% solve (Z'*Z - uv') x = b (assuming Z sparse, full column rank)

    zy = Z\(Z'\[u,b]);
    z = zy(:,1);
    y = zy(:,2);
    x = y + (v'*y)/(1-v'*z)*z;
    
return



function opts = setOpts(opts)
% This method sets the empty algorithmic options to their default values.
% --------------------------------------------------------------------------
% USAGE of "setOpts"
% opts = setOpts(opts)
% --------------------------------------------------------------------------
% INPUT
% opts:     custom algorithmic options
%
% OUTPUT
% opts:     complete algorithmic options
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------
           
    if ~isfield(opts, 'predLineSearch'); opts.predLineSearch = 1; end
    if ~isfield(opts, 'maxCorrSteps'); opts.maxCorrSteps = 4; end
    if ~isfield(opts, 'corrCheck'); opts.corrCheck = 1; end
    if ~isfield(opts, 'optimTol'); opts.optimTol = 1e-06; end
    if opts.optimTol < eps, opts.optimTol = eps; end
    if ~isfield(opts, 'debug'); opts.debug = 0; end
    if ~isfield(opts, 'maxCorrLSIters'); opts.maxCorrLSIters = 8; end
    if ~isfield(opts, 'maxPredSmallSteps'); opts.maxPredSmallSteps = 8; end
    if ~isfield(opts, 'maxItRefineSteps'); opts.maxItRefineSteps = 0; end
    if ~isfield(opts, 'verbose'); opts.verbose = 1; end
    if ~isfield(opts, 'maxIter'); opts.maxIter = 10000; end

return

function algParams = setAlgParams(gH_Params, opts)
% This method sets the algorithmic parameters.
% --------------------------------------------------------------------------
% USAGE of "setAlgParams"
% algParams = setAlgParams(gH_Params, opts)
% --------------------------------------------------------------------------
% INPUT
% gH_Params:	parameters associated with the method gH
% opts:         algorithmic options
%
% OUTPUT
% algParams:                        algorithmic parameters
% - algParams.maxIter:              maximum number of iterations
% - algParams.optimTol:             optimization tolerance parameter
% - algParams.alphaCorr:            corrector step size
% - algParams.predLSMulti:          predictor line search step size
%                                   multiplier
% - algParams.corrLSMulti:          corrector line search step size
%                                   multiplier
% - algParams.itRefineThreshold:    iterative refinement success threshold
% - algParams.maxCorrSteps:         maximum number of corrector steps
% - algParams.beta:                 large neighborhood parameter
% - algParams.eta:                  small neighborhood parameter
% - algParams.alphaPredLS:          initial predictor step size with line
%                                   search
% - algParams.alphaPredFix:         fixed predictor step size
% - algParams.alphaPred:            initial predictor step size
% - algParams.alphaPredThreshold:   minimum predictor step size
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    algParams.maxIter           = opts.maxIter;
    algParams.optimTol          = opts.optimTol;
    
    algParams.alphaCorr         = 1.0;
    algParams.predLSMulti       = 0.7;
    algParams.corrLSMulti       = 0.5;
    algParams.itRefineThreshold = 0.1;
    
    % parameters are chosen to make sure that each predictor 
    % step takes the current iterate from the eta-neighborhood to the
    % beta-neighborhood and each corrector phase takes the current 
    % iterate from the beta-neighborhood to the eta-neighborhood.
    % extra corrector steps are allowed to mitigate the effects of
    % finite precision.
    algParams.maxCorrSteps      = 2*opts.maxCorrSteps;

    % precomputed safe parameters
    switch opts.maxCorrSteps
        case 1
            if gH_Params.bnu < 10
                algParams.beta       = 0.1810;
                algParams.eta        = 0.0733;
                cPredFix             = 0.0225;
            elseif gH_Params.bnu < 100
                algParams.beta       = 0.2054;
                algParams.eta        = 0.0806;
                cPredFix             = 0.0263;              
            else
                algParams.beta       = 0.2190;
                algParams.eta        = 0.0836;
                cPredFix             = 0.0288;
            end
        case 2
            if gH_Params.bnu < 10
                algParams.beta       = 0.2084;
                algParams.eta        = 0.0502;
                cPredFix             = 0.0328;
            elseif gH_Params.bnu < 100
                algParams.beta       = 0.2356;
                algParams.eta        = 0.0544;
                cPredFix             = 0.0380;                  
            else
                algParams.beta       = 0.2506;
                algParams.eta        = 0.0558;
                cPredFix             = 0.0411;
            end
        case 4
            if gH_Params.bnu < 10
                algParams.beta       = 0.2387;
                algParams.eta        = 0.0305;
                cPredFix             = 0.0429;
            elseif gH_Params.bnu < 100
                algParams.beta       = 0.2683;
                algParams.eta        = 0.0327;
                cPredFix             = 0.0489;                  
            else
                algParams.beta       = 0.2844;
                algParams.eta        = 0.0332;
                cPredFix             = 0.0525;
            end
        otherwise
            error('The maximum number of corrector steps can be 1, 2, or 4.');
    end

    kx = algParams.eta + sqrt(2*algParams.eta^2 + gH_Params.bnu);
    algParams.alphaPredFix  = cPredFix/kx;
    algParams.alphaPredLS   = min(100 * algParams.alphaPredFix, 0.9999);
    algParams.alphaPredThreshold = (algParams.predLSMulti^opts.maxPredSmallSteps)*algParams.alphaPredFix;
        
    if opts.predLineSearch == 0
        % fixed predictor step size
        algParams.alphaPred   = algParams.alphaPredFix;
    else
        % initial predictor step size with line search
        algParams.alphaPred   = algParams.alphaPredLS;
    end

return

function soln = initSoln(x0, probData, gH, gH_Params)
% This method creates the central primal-dual iterate corresponding to x0.
% --------------------------------------------------------------------------
% USAGE of "initSoln"
% soln = initSoln(x0, probData, gH, gH_Params)
% --------------------------------------------------------------------------
% INPUT
% x0:           initial primal iterate
% probData:     data for the conic optimization problem
% gH:           method for computing the gradient and Hessian of the
%               barrier function
% gH_Params:    parameters associated with the method gH
%
% OUTPUT
% soln:         initial primal-dual iterate
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    [m, ~] = size(probData.A);
    [soln.in, soln.g, soln.H, soln.L] = gH(x0, gH_Params);
    soln.x       = x0;
    soln.y       = zeros(m, 1);
    soln.tau     = 1;
    soln.s       = -soln.g;
    soln.kappa   = 1;
    soln.mu      = (soln.x'*soln.s + soln.tau*soln.kappa)/gH_Params.bnu;

return

function [status, metrics] = term(soln, probData, algParams, termConsts)
% This method checks the termination criteria.
% --------------------------------------------------------------------------
% USAGE of "term"
% [status, statusString, metrics] = term(soln, probData, algParams, termConsts)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% algParams:    algorithmic parameters
% termConsts:   constants for termination criteria
%
% OUTPUT
% status:	    problem status code
% statusString: problem status as an interpretable string
% metrics:	    convergence metrics
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    x       = soln.x;
    tau     = soln.tau;
    s       = soln.s;
    kappa   = soln.kappa;
    y       = soln.y;
    mu      = soln.mu;
    
    A = probData.A;
    b = probData.b;
    c = probData.c;
    
    cx = c'*x;
    by = b'*y;
    
    % convergence metrics
    metrics.P = norm(A*x - tau*b,Inf) / termConsts.pRes;
    metrics.D = norm(A'*y + s - tau*c,Inf) / termConsts.dRes;
    metrics.G = abs(cx - by + kappa) / termConsts.comp;
    metrics.A = abs(cx - by) / (tau + abs(by));
    metrics.O = cx/tau;

    % complementarity gap of the initial iterate
    mu0 = 1;

    % termination criteria
    tol = algParams.optimTol;
    P   =    metrics.P   <= tol;
    D   =    metrics.D   <= tol;
    G   =    metrics.G   <= tol;
    AA  =    metrics.A   <= tol;
    T   =    tau <= tol * 1e-02 * max(1, kappa);
    K   =    tau <= tol * 1e-02 * min(1, kappa);
    M   =    mu  <= tol * 1e-02 * mu0;

    % are we at least in the ballpark?
    Papx   = metrics.P <= sqrt(tol);
    Dapx   = metrics.D <= sqrt(tol);
    Aapx  =  metrics.A <= sqrt(tol);
    
    if P && D && AA
        status = 1;         % success!
    elseif P && D && G && T
        if by > -tol && cx > -tol
            status = -1;    % P is infeasible
        elseif by < tol && cx < tol
            status = -2;    % D is infeasible
        elseif by > -tol && cx < tol
            status = -3;    % P and D both infeasible
        else
            status = -4;    % P or D is near infeasible
        end
    elseif K && M
        status = -8;        % ill-posed, likely no strong duality
    elseif Papx && Dapx && Aapx
        status = -6;        % no progress, but we are in the ballpark
    else
        status = -99;       % unknown error
    end

return

function results = prepResults(results, status, soln, probData, iter, time, options)
% This method prepares the final solution and iteration statistics.
% --------------------------------------------------------------------------
% USAGE of "prepResults"
% results = prepResults(results, soln, probData, iter)
% --------------------------------------------------------------------------
% INPUT
% results:              iteration statistics
% - results.alphaPred:  predictor step size at each iteration
% - results.betaPred:   neighborhood parameter at the end of the predictor
%                       phase at each iteration
% - results.etaCorr:    neighborhood parameter at the end of the corrector
%                       phase at each iteration
% - results.mu:         complementarity gap at each iteration
% status:               final problem status
% soln:                 current iterate
% probData:             data for the conic optimization problem
% iter:                 iteration count
% time:                 solver time
% options:              solver options structure
%
% OUTPUT
% results:	final solution and iteration statistics
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    results.nIterations     = iter;
    results.time            = time;
    
    % truncates arrays for iteration statistics
    results.alphaPred       = results.alphaPred(1:iter);
    results.betaPred        = results.betaPred(1:iter);
    results.etaCorr         = results.etaCorr(1:iter);
    results.mu              = results.mu(1:iter);
    
    % final solution
    results.x = soln.x/soln.tau;
    results.s = soln.s/soln.tau;
    results.y = soln.y/soln.tau;
    results.tau     = soln.tau;
    results.kappa   = soln.kappa;
    
    % final primal and dual objective values
    results.pObj    = probData.c'*results.x;
    results.dObj    = probData.b'*results.y;
    
    % solver options structure
    results.options = options;
    
    % final duality and complementarity gaps
    results.dGap    = results.pObj - results.dObj;
    results.cGap    = results.s'*results.x;
    
    % final relative duality and complementarity gaps
    results.rel_dGap = results.dGap/(1 + abs(results.pObj) + abs(results.dObj));
    results.rel_cGap = results.cGap/(1 + abs(results.pObj) + abs(results.dObj));
    
    % final primal and dual residuals
    results.pRes    = probData.b - probData.A*results.x;
    results.dRes    = probData.c - probData.A'*results.y - results.s;
    
    % final primal and dual infeasibilities
    results.pIn     = norm(results.pRes,Inf);
    results.dIn     = norm(results.dRes,Inf);
        
    % final relative primal and dual infeasibilities
    results.rel_pIn = results.pIn/(1+norm(probData.b,Inf));
    results.rel_dIn = results.dIn/(1+norm(probData.c,Inf));

    results.status = status;
    
    switch status
        case 1
            results.statusString = 'Optimal solution found.';
        case 0
            results.statusString = 'Number of iterations exceeded opts.maxIter.';
        case -1
            results.statusString = 'Primal infeasibility detected.';
        case -2
            results.statusString = 'Dual infeasibility detected.';
        case -3
            results.statusString = 'Primal and dual infeasibility detected.';
        case -4
            results.statusString = 'Problem is nearly primal or dual infeasible.';
        case -6
            results.statusString = 'Approximate/inaccurate solution found.';
        case -8
            results.statusString = 'Problem is ill-posed.';
        otherwise
            results.statusString = 'UNKNOWN';
    end

return

function [] = inputCheck(probData)
% This method checks the problem data for consistency.
% --------------------------------------------------------------------------
% USAGE of "inputCheck"
% inputCheck(probData)
% --------------------------------------------------------------------------
% INPUT
% probData: data for the conic optimization problem
%
% OUTPUT
% None.
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    [m, n] = size(probData.A);
    
    if m <= 0 || n <= 0
        error('Input matrix A must be nontrivial.');
    end
    if m ~= size(probData.b,1)
        error('Dimension of (column) vector b must match the number of rows in A.');
    end
    if n ~= size(probData.c,1)
        error('Dimension of (column) vector c must match the number of columns in A.');
    end
        
return

function say_hello(opts)

    if opts.verbose
        fprintf('\n*** alfonso (ver. 2020/07/20) by David Papp and Sercan Yildiz, (c) 2018.\n');
        if opts.predLineSearch
            fprintf('step size: line search,  ');
        else
            fprintf('step size: safe fixed,  ');
        end
        fprintf('optimality tolerance: %#.2e\n', opts.optimTol);
    end
    
return
