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
% Date: 
% 
% This code has been developed and tested with Matlab R2016b.
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
%                               Hessian of the barrier function
% gH_Params:                    parameters associated with the method gH
% - gH_Params.bnu:              complexity parameter of the augmented
%                               barrier (nu-bar)
% opts:                         algorithmic options
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
%                               value: 1e-06.
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
% --------------------------------------------------------------------------

    if nargin == 4
        opts = struct([]);
    elseif nargin < 4
        error('alfonso needs more input arguments.')
    elseif nargin > 5
        error('alfonso needs fewer input arguments.')
    end
    
    % sets algorithmic options
    opts = setOpts(opts); 
   
    % checks the problem data for consistency
    inputCheck(probData);
    
    [m, n] = size(probData.A);
    A = probData.A;
    b = probData.b;
    c = probData.c;
    
    probData.LHS = ...
    [ sparse(m,m)   A               -b            sparse(m,n)    sparse(m,1) ;
     -A'            sparse(n,n)      c           -speye(n)       sparse(n,1) ;
      b'           -c'               0            sparse(1,n)   -1          ;
      sparse(n,m)   speye(n)         sparse(n,1)  speye(n)       sparse(n,1) ;
      sparse(1,m)   sparse(1,n)      1            sparse(1,n)    1          ];

    % sets the solution method for the Newton system
    myLinSolve = @linSolve3;
    
    % sets algorithmic parameters
    algParams = setAlgParams(gH_Params, opts);
    
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
    
    termFlag = 0;
    numIters = 0;
    for iter = 1:algParams.maxIter

        % checks progress towards termination criteria
        [status, metrics] = term(soln, probData, algParams, termConsts);
                
        if termFlag || ~strcmpi(status, 'UNKNOWN')
            numIters = iter-1;
            disp(['Done in ', int2str(numIters), ' iterations.']);
            disp(['Status = ', status]);
            break;
        end
        
        % prints progress metrics
        if mod(iter,1)==0 && opts.verbose
            fprintf('%d: pObj=%d pIn=%d dIn=%d gap=%d tau=%d kap=%d mu=%d\n',...
            iter, metrics.O, metrics.P, metrics.D, metrics.A, soln.tau, soln.kappa, soln.mu);
        end

        % PREDICTOR PHASE
        [soln, alphaPred, betaPred, algParams, predStatus] =...
            pred(soln, probData, gH, gH_Params, myLinSolve, algParams, opts);
        results.alphaPred(iter) = alphaPred;
        results.betaPred(iter)  = betaPred;
        % raises a termination flag if predictor phase was not successful
        if predStatus == 0
            results.betaPred(iter) = results.etaCorr(iter-1);
            fprintf('Predictor could not improve the solution.\n');
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
                    fprintf('Corrector could not improve the solution.\n');
                    termFlag = 1; break;
                end
                % exits corrector phase if corrCheck == 1 and current
                % iterate is in the eta-neighborhood
                if (opts.corrCheck && corrIter < algParams.maxCorrSteps) || corrIter == algParams.maxCorrSteps
                    results.etaCorr(iter) = sqrt(sum((soln.L\soln.psi(1:end-1)).^2) +...
                        (soln.tau*soln.psi(end))^2)/soln.mu;
                    if results.etaCorr(iter) <= algParams.eta; break; end;
                end
            end
            % raises a termination flag if corrector phase was not successful
            if results.etaCorr(iter) > algParams.eta
                fprintf('Corrector phase finished outside the eta-neighborhood.\n');
                termFlag = 1;      
            end
        end
        results.mu(iter) = soln.mu;
    end
    
    % prepares final solution and iteration statistics
    results = prepResults(results, soln, probData, numIters);

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
                if opts.predLineSearch == 1; algParams.alphaPred = alpha; end;
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
                if opts.predLineSearch == 1; algParams.alphaPred = alpha; end;
                break;
            end
            % last two iterates were outside the beta-neighborhood and
            % alpha is very small
            if alpha < algParams.alphaPredThreshold
                predStatus  = 0; % predictor has failed
                alpha       = 0;
                betaAlpha   = Inf; % overwritten later in alfonso
                solnAlpha   = soln;
                if opts.predLineSearch == 1; algParams.alphaPred = alpha; end;
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
            LHS(end,m+n+1)              = soln.mu/soln.tau^2;
            
            exitFlag    = 0;
            res         = residual3p(LHS, delta, RHS);
            resNorm     = norm(res);  
            for iter = 1:opts.maxItRefineSteps
                if exitFlag; break; end;
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

function [delta, probData] = linSolve3(soln, probData, RHS)
% This method implements a block solver for the Newton system.
% --------------------------------------------------------------------------
% USAGE of "linSolve3"
% delta = linSolve3(soln, probData, RHS)
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
    
    ry     = RHS(1:m);
    rx     = RHS(m+(1:n));
    rtau   = RHS(m+n+1);
    rs     = RHS(m+n+1+(1:n));
    rkappa = RHS(end);
    
    Hic     = soln.L'\(soln.L\c);
    HiAt    = -soln.L'\(soln.L\A');
    Hirxrs  = soln.L'\(soln.L\(rx+rs));
    
    LHSdydtau   = [zeros(m), -b; b', soln.mu/soln.tau^2] - [A; -c']*[HiAt, Hic]/soln.mu;
    RHSdydtau   = [ry; rtau+rkappa] - [A; -c']*Hirxrs/soln.mu;
    dydtau      = LHSdydtau\RHSdydtau;
    dx          = (Hirxrs - [HiAt, Hic]*dydtau)/soln.mu;

    delta               = zeros(m+2*n+2, 1);
    delta(1:m)          = dydtau(1:m);
    delta(m+n+1)        = dydtau(m+1);
    delta(m+(1:n))      = dx;
    delta(m+n+1+(1:n))  = -rx - [A', -c]*dydtau;
    delta(end)          = -rtau + b'*dydtau(1:m) - c'*dx;
    
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
           
    if ~isfield(opts, 'predLineSearch'); opts.predLineSearch = 1; end;
    if ~isfield(opts, 'maxCorrSteps'); opts.maxCorrSteps = 4; end;
    if ~isfield(opts, 'corrCheck'); opts.corrCheck = 1; end;
    if ~isfield(opts, 'optimTol'); opts.optimTol = 1e-06; end;
    if ~isfield(opts, 'maxCorrLSIters'); opts.maxCorrLSIters = 8; end;
    if ~isfield(opts, 'maxPredSmallSteps'); opts.maxPredSmallSteps = 8; end;
    if ~isfield(opts, 'maxItRefineSteps'); opts.maxItRefineSteps = 0; end;    
    if ~isfield(opts, 'verbose'); opts.verbose = 1; end;
    
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
% - algParams.optimTol:               optimization tolerance parameter
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

    algParams.maxIter           = 10000;
    algParams.optimTol            = opts.optimTol;
    
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
% [status, metrics] = term(soln, probData, algParams, termConsts)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% algParams:    algorithmic parameters
% termConsts:   constants for termination criteria
%
% OUTPUT
% status:	problem status
% metrics:	convergence metrics
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
    P =    metrics.P   <= algParams.optimTol;
    D =    metrics.D   <= algParams.optimTol;
    G =    metrics.G   <= algParams.optimTol;
    AA =   metrics.A   <= algParams.optimTol;
    T =    tau <= algParams.optimTol * 1e-02 * max(1, kappa);
    K =    tau <= algParams.optimTol * 1e-02 * min(1, kappa);
    M =    mu  <= algParams.optimTol * 1e-02 * mu0;

    if P && D && AA
        status = 'Feasible and approximate optimal solution found.';
    elseif P && D && G && T
        status = 'Problem nearly primal or dual infeasible.';
    elseif K && M
        status = 'Problem is ill-posed.';
    else
        status = 'UNKNOWN';
    end

return

function results = prepResults(results, soln, probData, iter)
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
% soln:                 current iterate
% probData:             data for the conic optimization problem
% iter:                 iteration count
%
% OUTPUT
% results:	final solution and iteration statistics
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    results.nIterations     = iter;
    
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
    results.pIn     = norm(results.pRes);
    results.dIn     = norm(results.dRes);
        
    % final relative primal and dual infeasibilities
    results.rel_pIn = results.pIn/(1+norm(probData.b,Inf));
    results.rel_dIn = results.dIn/(1+norm(probData.c,Inf));
    
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
    if m ~= length(probData.b)
        error('Size of vector b must match the number of rows in A.');
    end
    if n ~= length(probData.c)
        error('Size of vector c must match the number of columns in A.');
    end
        
return
