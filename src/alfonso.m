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
%          Sercan Yildiz
%
% Version: 2024/07/15
% 
% This code has been developed and tested with Matlab R2023b.
% -------------------------------------------------------------------------

function results = alfonso(probData, x0, gH, gH_Params, opts)
% ALgorithm FOr Non-Symmetric Optimization
% This is the main method for the algorithm.
% --------------------------------------------------------------------------
% USAGE of "alfonso"
% results = alfonso(probData, x0, gH, gH_Params, opts)
% --------------------------------------------------------------------------
% INPUT
% probData:                     Data for the conic optimization problem.
% - probData.A:                 Constraint matrix.
% - probData.b:                 Right-hand side vector.
% - probData.c:                 Cost vector.
% x0:                           Initial primal iterate.
% gH:                           Method for computing the gradient and 
%                               the bilinear map corresponding to the
%                               Hessian inverse. (Function handle.)
% gH_Params:                    Parameters passed to gH; may be replaced
%                               by [] if no parameters needed.
%
% opts:                         Algorithmic options [default value]:
% - opts.maxIter:               Maximum number of interior-point iterations.
%                               [10000].
% - opts.linSolveFun:           Function handle for the Newton solver to use.
%                               (You may add your own!) [@linsolveA]
% - opts.predLineSearch:        0: fixed step size.
%                               1: line search.
%                               2: backtracking line search. [1]
% - opts.unsafeMult:            A multiplier for the safe fixed step size, 
%                               for those who like to live dangerously.
%                               [1]
% - opts.maxCorrSteps:          Maximum number of corrector steps. 
%                               Possible values: 1, 2, or 4. [4]
% - opts.corrCheck:             0: maxCorrSteps corrector steps are to be 
%                               performed at each corrector phase.
%                               1: corrector phase can be terminated before 
%                               maxCorrSteps corrector steps if the iterate 
%                               is in the eta-neighborhood. [1]
% - opts.optimTol:              Optimization tolerance parameter. [1e-06]
% - opts.maxCorrLSIters:        Maximum number of line search iterations in
%                               each corrector step. [8]
% - opts.maxSmallPredSteps:     Maximum number of predictor step size 
%                               reductions allowed with respect to the safe
%                               fixed step size. [8]
% - opts.verbose:               0: suppress output.
%                               1: progress is printed after each iteration. [1]
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
    sayHello(opts);
    
    % checks the problem data for consistency
    inputCheck(probData);

    stopwatch = tic;

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
    
    
    [~, n] = size(probData.A);
    A  = probData.A;
    b  = probData.b;
    c  = probData.c;
    At = A';
    probData.At = At;
    
    % if opts.maxItRefineSteps > 0 
    %     probData.LHS = ...
    %     [ sparse(m,m)   A               -b            sparse(m,n)    sparse(m,1) ;
    %      -A'            sparse(n,n)      c           -speye(n)       sparse(n,1) ;
    %       b'           -c'               0            sparse(1,n)   -1          ;
    %       sparse(n,m)   speye(n)         sparse(n,1)  speye(n)       sparse(n,1) ;
    %       sparse(1,m)   sparse(1,n)      1            sparse(1,n)    1          ];
    % end
    
    % sets the solution method for the Newton system
    % sets the solution method for the Newton system
    myLinSolve = opts.linSolveFun;
   
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
    termConsts.dRes = max([1, norm([At,speye(n),-c],Inf)]);
    termConsts.comp = max([1, norm([-c',b',1],Inf)]);

    % creates the central primal-dual iterate corresponding to x0
    soln = initSoln(x0, probData, gH, gH_Params);
    
    if ~opts.debug
        warning('off','MATLAB:nearlySingularMatrix');
    end
    termFlag  = 0;
    numIters  = 0;
    elapsed   = toc(stopwatch);
  
    for iter = 1:algParams.maxIter+1

        % checks progress towards termination criteria
        [status, metrics] = term(soln, probData, algParams, termConsts);
            
        % prints progress metrics
        if mod(iter,1)==0 && opts.verbose
            fprintf('%3d: pObj=%.6e pIn=%#.2e dIn=%#.2e gap=%#.2e tau=%#.2e kap=%#.2e mu=%.2e t=%#.2f s\n',...
                iter, metrics.O, metrics.P, metrics.D, metrics.A, soln.tau, soln.kappa, soln.mu, elapsed);
        end

        if termFlag || iter == algParams.maxIter+1 || (status ~= -99 && status ~= -6) 
            if iter == algParams.maxIter+1
                status = 0; % 'Number of iterations exceeded opts.maxIter.';
            end
            
            numIters = iter;
            break;
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
        else
            %?%soln.AHiAti = computeAHiAti(soln.Li, probData);
        end

        % CORRECTOR PHASE
        results.etaCorr(iter) = results.betaPred(iter); % Will be overwritten if corr phase is successful and it is measured.
        % Skips corrector phase if
        %   (corrCheck == 1 AND current iterate is already in the eta-neighborhood) OR
        %   termination flag was already raised
        if (~opts.corrCheck || results.etaCorr(iter) >= algParams.eta) && ~termFlag  % margin?
            for corrIter = 1:algParams.maxCorrSteps

                % A single corrector step, possibly with failsafe backtracking line search.
                % Fails only if could not stay inside the cone at all. (Does not check neighborhoods.)
                [soln, corrStatus] = corr(soln, probData, gH, gH_Params, myLinSolve, algParams, opts);
                %?%soln.AHiAti = computeAHiAti(soln.Li, probData);

                % Exits corrector phase and raises a termination flag if 
                % last corrector step was not successful. In this case, the
                % corrector ended outside the cone (not just the neighborhood).
                if corrStatus == 0
                    if opts.verbose
                        fprintf('Corrector could not improve the solution.\n');
                    end
                    termFlag = 1;
                    break;
                end
                % exits corrector phase if corrCheck == 1 and current
                % iterate is in the eta-neighborhood
                if (opts.corrCheck && corrIter < algParams.maxCorrSteps) || corrIter == algParams.maxCorrSteps
                    
                    if ~iscell(soln.Hi)
                        v = soln.psi(1:end-1);
                        psiHiPsi = v.'*soln.Hi(v);
                    else
                        idx = 0;
                        Kdims = probData.Kdims;
                        psiHiPsi = 0;
                        for k=1:length(Kdims)
                            v = soln.psi(idx+1:idx+Kdims(k));
                            psiHiPsi = psiHiPsi + v.'*soln.Hi{k}(v);
                            idx = idx+Kdims(k);
                        end
                    end
                    results.etaCorr(iter) = sqrt(psiHiPsi + (soln.tau*soln.psi(end))^2)/soln.mu;

                    if results.etaCorr(iter) < algParams.eta % ? margin
                        break;
                    end 

                end
            end

            % Raises a termination flag if corrector phase was not successful.
            % Unsafe fixed step size is exempted; in that case, all bets are off.
            if opts.debug && ~(opts.predLineSearch==0 && opts.unsafeMult>1) && results.etaCorr(iter) > algParams.eta
                if opts.verbose
                    fprintf('Corrector phase finished outside the eta-neighborhood.\n');
                end
                termFlag = 1;
            end
        end
        results.mu(iter) = soln.mu;
        
        elapsed = toc(stopwatch);
    end
    
    warning('on','MATLAB:nearlySingularMatrix');

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

    prevSolnAlpha = soln;  % At least we can return the curent point if everything else fails.

    x       = soln.x;
    tau     = soln.tau;
    s       = soln.s;
    kappa   = soln.kappa;
    y       = soln.y;
    
    A  = probData.A;
    At = probData.At;
    b  = probData.b;
    c  = probData.c;
    [m, n] = size(A);
    
    RHS = zeros(m+2*n+2,1);
    RHS(1:m+n+1)        = -[A*x - b*tau; -At*y + c*tau - s; b'*y - c'*x - kappa];
    RHS(m+n+2:end)      = -[s; kappa];
    
    % Computes the Newton direction.
    newtonDir = linSolveMain(soln, probData, RHS, myLinSolve);

    % Cannot even compute a predictor direction (severe numerical issue)
    if any(isnan(newtonDir.y))
        predStatus = 0;
        solnAlpha  = soln;
        alpha      = 0;
        betaAlpha  = NaN; % Not worth computing.
        return;
    end

    switch opts.predLineSearch
        case 0  % Fixed step size, no neighborhood checks.
            alpha            = algParams.alphaPred;

            solnAlpha.y      = y     + alpha*newtonDir.y;
            solnAlpha.x      = x     + alpha*newtonDir.x;
            solnAlpha.tau    = tau   + alpha*newtonDir.tau;
            solnAlpha.s      = s     + alpha*newtonDir.s;
            solnAlpha.kappa  = kappa + alpha*newtonDir.kappa;

            [solnAlpha.in, solnAlpha.g, solnAlpha.Hi, solnAlpha.Li] = gH(solnAlpha.x, gH_Params);

            if solnAlpha.in
                % Primal iterate is inside the cone, as it automatically should be.
                predStatus = 1;

                solnAlpha.mu	= (solnAlpha.x'*solnAlpha.s + solnAlpha.tau*solnAlpha.kappa)/gH_Params.bnu;
                solnAlpha.psi	= [solnAlpha.s; solnAlpha.kappa] + solnAlpha.mu*[solnAlpha.g; -1/solnAlpha.tau];
            else
                % This should not happen, except perhaps for numerical reasons.
                predStatus   = 0;

                % Restore the previous iterate.
                solnAlpha    = soln;
            end
            % Returned but not used outside of line search.
            betaAlpha    = NaN;

            return;

        case 1
            nSteps     = 0;
            alpha      = algParams.alphaPred;  % the initial step size or the last working one
            direction  = 0;  % +1 if increasing, -1 if decreasing, 0 if just started.

            % Loop until we find an acceptable step size or fail.
            while true
                nSteps = nSteps + 1;

                solnAlpha.y         = y     + alpha*newtonDir.y;
                solnAlpha.x         = x     + alpha*newtonDir.x;
                solnAlpha.tau       = tau   + alpha*newtonDir.tau;
                solnAlpha.s         = s     + alpha*newtonDir.s;
                solnAlpha.kappa     = kappa + alpha*newtonDir.kappa;

                [solnAlpha.in, solnAlpha.g, solnAlpha.Hi, solnAlpha.Li] = gH(solnAlpha.x, gH_Params);

                % Check if we are in the cone, and if so, if we are in the
                % beta neighborhood.
                if solnAlpha.in
                    % Primal iterate is inside the cone.
                   
                    % Compute the neighborhood radius.
                    solnAlpha.mu	= (solnAlpha.x'*solnAlpha.s +...
                        solnAlpha.tau*solnAlpha.kappa)/gH_Params.bnu;
                    solnAlpha.psi	= [solnAlpha.s; solnAlpha.kappa] +...
                        solnAlpha.mu*[solnAlpha.g; -1/solnAlpha.tau];

                    if ~iscell(solnAlpha.Hi)
                        v = solnAlpha.psi(1:end-1);
                        psiHiPsi = v.'*solnAlpha.Hi(v);
                    else
                        idx = 0;
                        Kdims = probData.Kdims;
                        psiHiPsi = 0;
                        for k=1:length(Kdims)
                            v = solnAlpha.psi(idx+1:idx+Kdims(k));
                            psiHiPsi = psiHiPsi + v.'*solnAlpha.Hi{k}(v);
                            idx = idx+Kdims(k);
                        end
                    end
                    % The neighborhood radius beta after the alpha step:
                    betaAlpha = sqrt(psiHiPsi + (solnAlpha.tau*solnAlpha.psi(end))^2)/solnAlpha.mu;
                    if (betaAlpha < algParams.beta)  % ? margin
                        inNbhood = true;
                    else
                        inNbhood = false;
                    end
                else
                    % Not in the neighborhood
                    inNbhood = false;
                    %solnAlpha.mu = NaN;
                end

                if inNbhood
                    % Good step size.

                    % Save the last good iterate.
                    prevSolnAlpha = solnAlpha;
                    prevBetaAlpha = betaAlpha;
                    prevAlpha     = alpha;

                    % Reached maximum step size?
                    maxStepSize = 0.99;
                    if alpha >= maxStepSize
                        %solnAlpha.mu    = (solnAlpha.x'*solnAlpha.s + solnAlpha.tau*solnAlpha.kappa)/gH_Params.bnu;
                        break;
                    end

                    if direction >= 0
                        % Increase the step size if possible.
                        direction = 1;  % If we weren't increasing yet, we are now.

                        % Can we increase?
                        if alpha / algParams.predLS1Multi < maxStepSize
                            alpha = alpha / algParams.predLS1Multi;
                            continue;
                        else
                            alpha = maxStepSize;
                            continue;
                        end
                    else
                        % We have been decreasing, but found a good step size.
                        break;
                    end
                else
                    % Step size too large.
                    if direction == 1
                        % We have been increasing, but we went too far.
                        % Accept previous stepsize.
                        solnAlpha = prevSolnAlpha;
                        betaAlpha = prevBetaAlpha;
                        alpha     = prevAlpha;
                        break;

                    else
                        direction = -1;  % If we weren't decreasing yet, we are now.

                        % Lower step size if possible.
                        alpha = alpha * algParams.predLS1Multi;
                        % Have we failed?
                        if alpha < algParams.alphaPredThreshold
                            predStatus  = 0; % predictor has failed
                            solnAlpha = prevSolnAlpha;
                            betaAlpha = alpha;
                            return;
                        end
                    end
                end

            end % end while
            
            % If we are here, we found the good step size.
            predStatus = 1;
            % Will start from here next time.
            algParams.alphaPred = alpha;

            return;

        case 2 %
            nSteps     = 0;
            alpha      = algParams.alphaPred;  % the initial step size or the last working one

            % Loop until we find an acceptable step size or fail.
            while true
                nSteps = nSteps + 1;

                solnAlpha.y         = y     + alpha*newtonDir.y;
                solnAlpha.x         = x     + alpha*newtonDir.x;
                solnAlpha.tau       = tau   + alpha*newtonDir.tau;
                solnAlpha.s         = s     + alpha*newtonDir.s;
                solnAlpha.kappa     = kappa + alpha*newtonDir.kappa;

                [solnAlpha.in, solnAlpha.g, solnAlpha.Hi, solnAlpha.Li] = gH(solnAlpha.x, gH_Params);

                if ~solnAlpha.in
                    % Primal iterate is outside the cone.
                    % Reduce step size and repeat if possible.

                    alpha = alpha * algParams.predLS2Multi;
                    % Have we failed?
                    if alpha < algParams.alphaPredThreshold
                        predStatus  = 0; % predictor has failed
                        %alpha       = 0;
                        % Restore the previous iterate.
                        solnAlpha    = soln;
                        % Reported but not used.
                        betaAlpha    = NaN;
                        return;
                    end

                else
                    % Compute the neighborhood radius.

                    solnAlpha.mu    = (solnAlpha.x'*solnAlpha.s + solnAlpha.tau*solnAlpha.kappa)/gH_Params.bnu;
                    solnAlpha.psi	= [solnAlpha.s; solnAlpha.kappa] + solnAlpha.mu*[solnAlpha.g; -1/solnAlpha.tau];

                    if ~iscell(solnAlpha.Hi)
                        v = solnAlpha.psi(1:end-1);
                        psiHiPsi = v.'*solnAlpha.Hi(v);
                    else
                        idx = 0;
                        Kdims = probData.Kdims;
                        psiHiPsi = 0;
                        for k=1:length(Kdims)
                            v = solnAlpha.psi(idx+1:idx+Kdims(k));
                            psiHiPsi = psiHiPsi + v.'*solnAlpha.Hi{k}(v);
                            idx = idx+Kdims(k);
                        end
                    end
                    % The neighborhood radius beta after the alpha step:
                    betaAlpha = sqrt(psiHiPsi + (solnAlpha.tau*solnAlpha.psi(end))^2)/solnAlpha.mu;

                    % If we are in the beta-neighborhood, done.
                    % Otherwise reduce step size and continue if possible.
                    if (betaAlpha < algParams.beta)
                        predStatus = 1;
                        %algParams.alphaPred = alpha;  % This might need more testing to see if it's worth it. ??
                        return;
                    else
                        alpha = alpha * algParams.predLS2Multi;
                        % Have we failed?
                        if alpha < algParams.alphaPredThreshold
                            predStatus  = 0; % predictor has failed
                            %alpha       = 0;
                            return;
                        end
                    end
                end % end if ~solnAlpha.in
            end  % endwhile

    end % end switch
        
return

function [solnAlpha, corrStatus] = corr(soln, probData, gH, gH_Params, myLinSolve, algParams, opts)
% This method performs a single corrector step in the algorithm.
% It does not do neighborhood checks, which is included in the main script.
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
    
    dsoln = linSolveMain(soln, probData, RHS, myLinSolve);

    % Cannot even compute a corrector direction (severe numerical issue).
    if any(isnan(dsoln.y))
        corrStatus = 0;
        solnAlpha = soln;
        return;
    end
    
    alpha = algParams.alphaCorr;
    corrStatus = 1;
    
    % does line search to make sure next primal iterate remains inside the cone
    for nSteps = 1:opts.maxCorrLSIters
        
        solnAlpha.y     = soln.y     + alpha*dsoln.y;
        solnAlpha.x     = soln.x     + alpha*dsoln.x;
        solnAlpha.tau   = soln.tau   + alpha*dsoln.tau;
        solnAlpha.s     = soln.s     + alpha*dsoln.s;
        solnAlpha.kappa = soln.kappa + alpha*dsoln.kappa;
        
        [solnAlpha.in, solnAlpha.g, solnAlpha.Hi, solnAlpha.Li] = gH(solnAlpha.x, gH_Params);
        
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

function dsoln = linSolveMain(soln, probData, RHS, myLinSolve)
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
%
% OUTPUT
% dsoln:	computed Newton direction
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------
    
    [m, n] = size(probData.A);
    delta  = myLinSolve(soln, probData, RHS);
   
    dsoln.y     = delta(1:m);
    dsoln.x     = delta(m+(1:n));
    dsoln.tau   = delta(m+n+1);
    dsoln.s     = delta(m+n+1+(1:n));
    dsoln.kappa = delta(end);
    
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
    if ~isfield(opts, 'unsafeMult'); opts.unsafeMult = 1.0; end
    if ~isfield(opts, 'maxCorrSteps'); opts.maxCorrSteps = 4; end
    if ~isfield(opts, 'corrCheck'); opts.corrCheck = 1; end
    if ~isfield(opts, 'optimTol'); opts.optimTol = 1e-06; end
    if opts.optimTol < eps, opts.optimTol = eps; end
    if ~isfield(opts, 'debug'); opts.debug = 0; end
    if ~isfield(opts, 'maxCorrLSIters'); opts.maxCorrLSIters = 8; end
    if ~isfield(opts, 'maxPredSmallSteps'); opts.maxPredSmallSteps = 8; end
    if ~isfield(opts, 'verbose'); opts.verbose = 1; end
    if ~isfield(opts, 'maxIter'); opts.maxIter = 10000; end
    if ~isfield(opts, 'linSolveFun'); opts.linSolveFun = @linsolveA; end

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
% - algParams.predLS1Multi:         predictor line search 1 step size
%                                   multiplier
% - algParams.predLS2Multi:         predictor line search 2 step size
%                                   multiplier
% - algParams.corrLSMulti:          corrector line search step size
%                                   multiplier
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
    algParams.predLS1Multi      = 0.7;
    algParams.predLS2Multi      = 0.7;
    algParams.corrLSMulti       = 0.5;
    
    % parameters are chosen to make sure that each predictor 
    % step takes the current iterate from the eta-neighborhood to the
    % beta-neighborhood and each corrector phase takes the current 
    % iterate from the beta-neighborhood to the eta-neighborhood.
    % extra corrector steps are allowed to mitigate the effects of
    % finite precision.
    algParams.maxCorrSteps      = 10*opts.maxCorrSteps;

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
    algParams.alphaPredLS   = min(100 * algParams.alphaPredFix, 0.99);
    
    switch opts.predLineSearch
        case 0
            % Fixed predictor step size.
            % The safe step size is usually "too safe".
            algParams.alphaPred   = opts.unsafeMult*algParams.alphaPredFix;
        case 1
            % Initial and minimum predictor step size for full line search.
            algParams.alphaPred   = algParams.alphaPredLS;
            algParams.alphaPredThreshold = 1e-04;
        case 2
            % Initial and minimum predictor step size for backtracking line search.
            algParams.alphaPred   = 0.99;
            algParams.alphaPredThreshold = 1e-04;
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
% soln:         initial primal-dual iteratesoln = str
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    [m, ~] = size(probData.A);
    [in, g, Hi, Li] = gH(x0, gH_Params);
    %?%AHiAti = computeAHiAti(Li, probData);

    % Packing it up for the soln structure.
    if iscell(Li)
        Li = {Li};
        Hi = {Hi};
    end

    soln = struct(          ...
        'in',  in,          ...
        'g',   g,           ...
        'Li',  Li,         ...
        'Hi',  Hi,          ...,
        'x',   x0,          ...
        'y',   zeros(m, 1), ...
        'tau', 1,           ...
        's',   -g,          ...
        'kappa', 1,         ...
        'mu', 1             ...
    );       

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
    
    A  = probData.A;
    At = probData.At;
    b  = probData.b;
    c  = probData.c;
    
    cx = c'*x;
    by = b'*y;
    
    % convergence metrics
    metrics.P = norm(A*x - tau*b,Inf) / termConsts.pRes;
    metrics.D = norm(At*y + s - tau*c,Inf) / termConsts.dRes;
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

function results = prepResults(results, status, soln, probData, iter, time, opts)
% This method prepares the final solution and iteration statistics.
% --------------------------------------------------------------------------
% USAGE of "prepResults"
% results = prepResults(results, soln, probData, iter, time, opts)
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
% statusString:         final problem status
% soln:                 current iterate
% probData:             data for the conic optimization problem
% iter:                 iteration count
% time:                 elapsed time in seconds
% opts:                 alfonso options
%
% OUTPUT
% results:	final solution and iteration statistics
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    results.nIterations     = iter;
    results.time            = time;
    
    % Truncates or removes the arrays for iteration statistics.
    if opts.debug
        results.alphaPred       = results.alphaPred(1:iter);
        results.betaPred        = results.betaPred(1:iter);
        results.etaCorr         = results.etaCorr(1:iter);
        results.mu              = results.mu(1:iter);
    else
        results = rmfield(results,{'alphaPred','betaPred','etaCorr','mu'});
    end
    
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
    
    % final primal and dual linear infeasibilities
    results.pRes    = probData.b - probData.A*results.x;
    results.dRes    = probData.c - probData.At*results.y - results.s;

    % final primal and dual infeasibilities
    results.pIn     = norm(results.pRes);
    results.dIn     = norm(results.dRes);
        
    % final relative primal and dual linear infeasibilities
    results.rel_pIn = results.pIn/(1+norm(probData.b,Inf));
    results.rel_dIn = results.dIn/(1+norm(probData.c,Inf));

    % final relative duality and complementarity gaps
    results.rel_dGap = results.dGap/(1 + abs(results.pObj) + abs(results.dObj));
    results.rel_cGap = results.cGap/(1 + abs(results.pObj) + abs(results.dObj));

    results.status = status;
    
    switch status
        case 0
            results.statusString = 'Number of iterations exceeded opts.maxIter.';
        case 1
            results.statusString = 'Optimal solution found.';
        case -1
            results.statusString = 'Primal infeasibility detected.';
        case -2
            results.statusString = 'Dual infeasibility detected.';
        case -3
            results.statusString = 'Primal and dual infeasibility detected.';
        case -4
            results.statusString = 'Problem is nearly primal or dual infeasible.';
        case -6
            results.statusString = 'Inaccurate solution found.';
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

function sayHello(opts)

    if opts.verbose
        fprintf('\n*** alfonso (ver. 2024/07/15) by David Papp and Sercan Yildiz, (c) 2018.\n');
        switch opts.predLineSearch
            case 0
                fprintf('step size: safe fixed,  ');
            case 1
                fprintf('step size: line search,  ');
            case 2
                fprintf('step size: backtracking,  ');
        end
        fprintf('optimality tolerance: %#.2e\n', opts.optimTol);
    end
    
return
