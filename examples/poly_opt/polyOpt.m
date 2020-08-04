% This code is an implementation of the sum-of-squares optimization approach 
% based on non-symmetric conic optimization and polynomial interpolants 
% presented in:
%
% D. Papp and S. Yildiz. Sum-of-squares optimization without semidefinite 
% programming. Available at https://arxiv.org/abs/1712.01792.
%
% The implementation formulates and solves the polynomial optimization
% problems described in the same reference. 
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
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% None.
% -------------------------------------------------------------------------

function sol = polyOpt(intParams, polyName, tol)
% This is the main method for the sum-of-squares optimization approach to 
% to the polynomial optimization problem.
% --------------------------------------------------------------------------
% USAGE of "polyOpt"
% results = polyOpt(intParams, polyName, tol)
% --------------------------------------------------------------------------
% INPUT
% intParams:        data for the interpolant basis (as generated in
%                   ChebInterval.m, for instance)
% - intParams.n:    number of arguments to the polynomials
% - intParams.d:    largest degree of polynomials to be squared
% - intParams.L:    dimension of the space of (intParams.n)-variate
%                   degree-(intParams.d) polynomials
% - intParams.U:    dimension of the space of (intParams.n)-variate
%                   degree-(2*intParams.d) polynomials
% - intParams.pts:  interpolation points. (intParams.U x intParams.n) array. 
% - intParams.w:    quadrature weights for the interpolation points
% - intParams.P0:   evaluations of a basis for the space of
%                   (intParams.n)-variate degree-(intParams.d) 
%                   polynomials at the interpolation points.
%                   (intParams.U x intParams.L) array. columns are indexed
%                   with the basis polynomials. it is assumed that the first 
%                   nchoosek(intParams.n+k,intParams.n) columns are a
%                   basis for the space of (intParams.n)-variate degree-k
%                   polynomials for k = 1,...,intParams.d.
% - intParams.P:    similar to intParams.P0
% polyName:         name of polynomial to be minimized. see the function
%                   setPolyParams below for options.
% tol:              tolerance parameter for optimization
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
    
    n   = intParams.n;
    d   = intParams.d;
    U   = intParams.U;
    L   = intParams.L;
    P0  = intParams.P0;
    pts = intParams.pts;

    [polyDeg, lb, ub] = setPolyParams(polyName, n);

    if d < ceil(polyDeg/2)
        error(strcat(polyName,' requires d >= ',num2str(ceil(polyDeg/2)),'.'));
    end
    
    % transforms points to fit the domain
    scale   = (ub-lb)/2;
    shift   = (lb+ub)/2;
    pts     = bsxfun(@plus,bsxfun(@times,pts,scale'),shift');    
    wtVals  = bsxfun(@minus,pts,lb').*bsxfun(@minus,ub',pts);

    % ORDER OF VARIABLES 
    % x \in WSOS_(n,2*d)^*
   
    % WEIGHTS
    % weight #0: 1
    % degree of associated SOS multiplier: 2*d
    % weight #j for j = 1,...,n: (lb_j+t_j)(ub_j-t_j)
    % degree of associated SOS multiplier: 2*d-2

    LWts = nchoosek(n+d-1,n)*ones(n,1);
    
    % PARAMETERS ASSOCIATED WITH THE METHOD gH_polyOpt
    % see input description for gH_polyOpt for details

    gH_Params.n = n;
    gH_Params.d = d;
    gH_Params.U = U;
    
    gH_Params.L     = L;
    gH_Params.LWts  = LWts;
    nu              = L + sum(LWts);
    gH_Params.bnu	= nu+1;

    % P0 has columns of Chebyshev polynomial evaluations:
    gH_Params.P = P0;
    % associated positive semidefinite cone constraint:
    % P0'*diag(x)*P0 >= 0
    PWts = cell(n,1);
    for j = 1:n
        PWts{j} = diag(sqrt(wtVals(:,j)))*P0(:,1:LWts(j));
        % associated positive semidefinite cone constraint: 
        % PWts{j}'*diag(x)*PWts{j} >= 0
    end
    gH_Params.PWts = PWts;
    
    % x \in WSOS_(n,2*d)^* <=> 
    % P'*diag(x)*P >= 0, PWts{1}'*diag(x)*PWts{1} >= 0, PWts{2}'*diag(x)*PWts{2} >= 0,...

    % DATA FOR THE CONIC OPTIMIZATION PROBLEM
    probData.A = ones(1,U);
    probData.b = 1;
    probData.c = setPolyVals(polyName, pts);

    tic
    
    % INITIAL PRIMAL ITERATE
    x0 = ones(U,1);
    [~, g0, ~, ~] = gH_polyOpt(x0, gH_Params);
    % scaling factor for the primal problem
    rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
    % scaling factor for the dual problem
    rD = max((1+abs(g0))./(1+abs(probData.c)));
    % initial primal iterate
    x0 = repmat(sqrt(rP*rD), U, 1);
    
    % CUSTOM ALGORITHMIC OPTIONS
    opts.optimTol = tol;
    
    % CALL TO alfonso
    sol = alfonso(probData,x0,@gH_polyOpt,gH_Params,opts);
    fprintf('alfonso is done.\n');
    
    toc
    
    % prints error statistics
    fprintf('\n');
    fprintf('FINAL:\n');
    fprintf('Relative primal infeasibility: %d\n', sol.rel_pIn);
    fprintf('Relative dual infeasibility: %d\n', sol.rel_dIn);
    fprintf('Relative duality gap: %d\n', sol.rel_dGap);
    fprintf('Relative complementarity gap: %d\n\n', sol.rel_cGap);
    
return

function [in, g, H, L] = gH_polyOpt(x, params)
% This method computes the gradient and Hessian of the barrier function for
% the polynomial optimization problem.
% --------------------------------------------------------------------------
% USAGE of "gH_polyOpt"
% [in, g, H, L] = gH_polyOpt(x, params)
% --------------------------------------------------------------------------
% INPUT
% x:                    primal iterate
% params:               parameters associated with the method gH
% - params.n:           number of arguments to the polynomials
% - params.d:           largest degree of polynomials to be squared
% - params.U:           dimension of the space of (params.n)-variate 
%                       degree-(2*params.d) polynomials
% - params.L:           dimension of the space of (params.n)-variate 
%                       degree-(params.d) polynomials.
% - params.LWts:        dimensions of the "weighted" polynomial spaces. 
%                       params.LWts(j) is the dimension of the space of
%                       (params.n)-variate degree-(params.d-1) polynomials
%                       for j = 1,...,n.
% - params.bnu:         complexity parameter of the augmented barrier (nu-bar)
% - params.P:           evaluations of the basis for the space of 
%                       (params.n)-variate degree-(params.d) polynomials
%                       at the interpolation points
% - params.PWts:        evaluations of "weighted" basis polynomials at the
%                       interpolation points. params.PWts{j} is
%                       the evaluations of the basis for the "weighted"
%                       space of (params.n)-variate degree-(params.d-1)
%                       polynomials at the interpolation points for 
%                       j = 1,...,n. the weight corresponding to 
%                       params.PWts{j} is sqrt((lb_j+t_j)(ub_j-t_j)).
%
% OUTPUT
% in:	0 if x is not in the interior of the cone. 1 if x is in the 
%       interior of the cone.
% g:    gradient of the barrier function at x
% H:    Hessian of the barrier function at x
% L:    Cholesky factor of H
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    n       = params.n;
    P       = params.P;
    PWts    = params.PWts;
    
    % ORDER OF VARIABLES 
    % x \in WSOS_(n,2*d)^*
    
    % for the weight 1
    [in, g, H] = gH_SOSWt(x,P);
    
    if in == 1
        for j = 1:n
            % for the weight (lb_j+t_j)(ub_j-t_j)
            [inWt, gWt, HWt] = gH_SOSWt(x,PWts{j});
            in  = in & inWt;
            if in == 1
                g   = g+gWt;
                H   = H+HWt;
            else
                g   = NaN;
                H   = NaN;
                break;
            end
        end
    end
    
    if in == 1        
        % checks positive semidefiniteness of H one more time.
        % H may fail Cholesky factorization even if in == 1
        % due to numerical errors in summing up HWt's above.
        [L, err] = chol(H,'lower');
        in = in & (err == 0);
    else
        L = NaN;
    end
    
return

function [in, g, H] = gH_SOSWt(x, P)
% This method computes the gradient and Hessian of a barrier function term
% corresponding to a single weight for the polynomial optimization problem.
% --------------------------------------------------------------------------
% USAGE of "gH_SOSWt"
% [in, g, H] = gH_SOSWt(x, P)
% --------------------------------------------------------------------------
% INPUT
% x:    primal iterate
% P:    evaluations of "weighted" basis polynomials at the
%       interpolation points
%
% OUTPUT
% in:   0 if P'*diag(x)*P is not positive definite. 1 if P'*diag(x)*P is
%       positive definite.
% g:    gradient of the barrier function term corresponding to a single
%       weight at x
% H:    Hessian of the barrier function term corresponding to a single
%       weight at x
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    Y = P'*diag(x)*P;

    if ~issymmetric(Y)
        Y = (Y+Y')/2;
    end

    [L,err] = chol(Y,'lower');
    if err > 0
        in = 0;
        g = NaN;
        H = NaN;
    else
        in = 1;
        V = L\P';
        VtV = V'*V;

        g = -diag(VtV);
        H = VtV.^2;
    end
    
return

function [polyDeg, lb, ub] = setPolyParams(polyName, n)
% This method sets the values of certain parameters associated with the
% polynomial optimization problem.
% --------------------------------------------------------------------------
% USAGE of "setPolyParams"
% [polyDeg, lb, ub] = setPolyParams(polyName, n)
% --------------------------------------------------------------------------
% INPUT
% polyName:	name of polynomial to be minimized
% n:        number of arguments to the polynomials
%
% OUTPUT
% polyDeg:	degree of polynomial to be minimized
% lb:       lower bounds defining the polynomial's standard domain
% ub:       upper bounds defining the polynomial's standard domain
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    switch polyName
        case 'robinson'            
            if n ~= 2; error('robinson requires 2 arguments.'); end
            polyDeg = 6;
            lb = repmat(-1, n, 1);
            ub = repmat(1, n, 1);
        case 'rosenbrock'
            if n ~= 2; error('rosenbrock requires 2 arguments.'); end
            polyDeg = 4;
            lb = repmat(-1, n, 1);
            ub = repmat(1, n, 1);
        case 'motzkin'
            if n ~= 2; error('motzkin requires 2 arguments.'); end
            polyDeg = 6;
            lb = repmat(-1, n, 1);
            ub = repmat(1, n, 1);
        case 'schwefel'
            if n ~= 3; error('schwefel requires 3 arguments.'); end
            polyDeg = 4;  
            lb = repmat(-10, n, 1);
            ub = repmat(10, n, 1);
        case 'reaction-diffusion'
            if n ~= 3; error('reaction-diffusion requires 3 arguments.'); end
            polyDeg = 2;  
            lb = repmat(-5, n, 1);
            ub = repmat(5, n, 1);
        case 'caprasse'    
            if n ~= 4; error('caprasse requires 4 arguments.'); end
            polyDeg = 4;            
            lb = repmat(-0.5, n, 1);
            ub = repmat(0.5, n, 1);
        case 'lotka-volterra'
            if n ~= 4; error('lotka-volterra requires 4 arguments.'); end
            polyDeg = 3;
            lb = repmat(-2, n, 1);
            ub = repmat(2, n, 1);
        case 'butcher'
            if n ~= 6; error('butcher requires 6 arguments.'); end
            polyDeg = 3;
            lb = [-1; -0.1; -0.1; -1; -0.1; -0.1];         
            ub = [0; 0.9; 0.5; -0.1; -0.05; -0.03];
        case 'magnetism7'
            if n ~= 7; error('magnetism7 requires 7 arguments.'); end
            polyDeg = 2;            
            lb = repmat(-1, n, 1);
            ub = repmat(1, n, 1);
        case 'heart'
            if n ~= 8; error('heart requires 8 arguments.'); end
            polyDeg = 4;
            lb = [-0.1; 0.4; -0.7; -0.7; 0.1; -0.1; -0.3; -1.1];         
            ub = [0.4; 1; -0.4; 0.4; 0.2; 0.2; 1.1; -0.3];
        otherwise
            error('The polynomial name is not recognized.');
    end
    %#ok<*REPMAT>

return

function vals = setPolyVals(polyName, pts)
% This method computes evaluations of the polynomial to be optimized at
% the interpolation points.
% --------------------------------------------------------------------------
% USAGE of "setPolyVals"
% vals = setPolyVals(polyName, pts)
% --------------------------------------------------------------------------
% INPUT
% polyName:	name of polynomial to be minimized
% pts:      interpolation points
%
% OUTPUT
% vals: evaluations of the polynomial to be optimized at
%       the interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    switch polyName
        case 'robinson'
            vals = robinsonVals(pts);
        case 'rosenbrock'
            vals = rosenbrockVals(pts);
        case 'motzkin'
            vals = motzkinVals(pts);
        case 'schwefel'
            vals = schwefelVals(pts);     
        case 'reaction-diffusion'
            vals = reactionDiffusionVals(pts);  
        case 'caprasse'
            vals = caprasseVals(pts);
        case 'lotka-volterra'
            vals = lotkaVolterraVals(pts);
        case 'butcher'
            vals = butcherVals(pts);
        case 'magnetism7'
            vals = magnetism7Vals(pts);            
        case 'heart'
            vals = heartVals(pts);
    end
    
return

function vals = robinsonVals(pts)
% This method computes evaluations of Robinson's polynomial at the
% interpolation points.
% --------------------------------------------------------------------------
% USAGE of "robinsonVals"
% vals = robinsonVals(pts)
% --------------------------------------------------------------------------
% INPUT
% pts:	interpolation points
%
% OUTPUT
% vals: evaluations of Robinson's polynomial at the interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    vals = 1 + pts(:,1).^6 + pts(:,2).^6 - (pts(:,1).^4).*(pts(:,2).^2) +...
        pts(:,1).^4 - (pts(:,2).^4).*(pts(:,1).^2) + pts(:,2).^4 -...
        pts(:,1).^2 + pts(:,2).^2 + 3*(pts(:,1).^2).*(pts(:,2).^2);

return

function vals = rosenbrockVals(pts)
% This method computes evaluations of Rosenbrock's polynomial at the
% interpolation points.
% --------------------------------------------------------------------------
% USAGE of "rosenbrockVals"
% vals = rosenbrockVals(pts)
% --------------------------------------------------------------------------
% INPUT
% pts:	interpolation points
%
% OUTPUT
% vals: evaluations of Rosenbrock's polynomial at the interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    vals = 1 - 2*pts(:,1) + pts(:,1).^2 + 100*pts(:,1).^4 -...
        200*(pts(:,1).^2).*pts(:,2) + 100*pts(:,2).^2;

return
    
function vals = motzkinVals(pts)
% This method computes evaluations of Motzkin's polynomial at the
% interpolation points.
% --------------------------------------------------------------------------
% USAGE of "motzkinVals"
% vals = motzkinVals(pts)
% --------------------------------------------------------------------------
% INPUT
% pts:	interpolation points
%
% OUTPUT
% vals: evaluations of Motzkin's polynomial at the interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    vals = 1 - 48*(pts(:,1).^2).*(pts(:,2).^2) + 64*(pts(:,1).^2).*(pts(:,2).^4) +...
        64*(pts(:,1).^4).*(pts(:,2).^2);

return

function vals = schwefelVals(pts)
% This method computes evaluations of Schwefel's polynomial at the
% interpolation points.
% --------------------------------------------------------------------------
% USAGE of "schwefelVals"
% vals = schwefelVals(pts)
% --------------------------------------------------------------------------
% INPUT
% pts:	interpolation points
%
% OUTPUT
% vals: evaluations of Schwefel's polynomial at the interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    vals = (pts(:,1) - pts(:,2).^2).^2 + (pts(:,2) - 1).^2 +...
        (pts(:,1) - pts(:,3).^2).^2 + (pts(:,3) - 1).^2;

return

function vals = reactionDiffusionVals(pts)
% This method computes evaluations of the 3-variable reaction-diffusion
% polynomial at the interpolation points.
% --------------------------------------------------------------------------
% USAGE of "reactionDiffusionVals"
% vals = reactionDiffusionVals(pts)
% --------------------------------------------------------------------------
% INPUT
% pts:	interpolation points
%
% OUTPUT
% vals: evaluations of the 3-variable reaction-diffusion polynomial at the 
%       interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    vals = -pts(:,1) + 2*pts(:,2) - pts(:,3) - 0.835634534*pts(:,2).*(1+pts(:,2));

return

function vals = caprasseVals(pts)
% This method computes evaluations of Caprasse's system at the interpolation 
% points.
% --------------------------------------------------------------------------
% USAGE of "caprasseVals"
% vals = caprasseVals(pts)
% --------------------------------------------------------------------------
% INPUT
% pts:	interpolation points
%
% OUTPUT
% vals: evaluations of Caprasse's system at the interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    vals = -pts(:,1).*(pts(:,3).^3) + 4*pts(:,2).*(pts(:,3).^2).*pts(:,4) +...
        4*pts(:,1).*pts(:,3).*(pts(:,4).^2) + 2*pts(:,2).*(pts(:,4).^3) +...
        4*pts(:,1).*pts(:,3) + 4*(pts(:,3).^2) - 10*pts(:,2).*pts(:,4) -...
        10*(pts(:,4).^2) + 2;

return

function vals = lotkaVolterraVals(pts)
% This method computes evaluations of the adaptive Lotka-Volterra system at 
% the interpolation points.
% --------------------------------------------------------------------------
% USAGE of "lotkaVolterraVals"
% vals = lotkaVolterraVals(pts)
% --------------------------------------------------------------------------
% INPUT
% pts:	interpolation points
%
% OUTPUT
% vals: evaluations of the adaptive Lotka-Volterra system at the 
%       interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    vals = pts(:,1).*(pts(:,2).^2) + pts(:,1).*(pts(:,3).^2) +...
        pts(:,1).*(pts(:,4).^2) - 1.1*pts(:,1) + 1;
    
return

function vals = butcherVals(pts)
% This method computes evaluations of Butcher's polynomial at the interpolation
% points.
% --------------------------------------------------------------------------
% USAGE of "butcherVals"
% vals = butcherVals(pts)
% --------------------------------------------------------------------------
% INPUT
% pts:	interpolation points
%
% OUTPUT
% vals: evaluations of Butcher's polynomial at the interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    vals = pts(:,6).*(pts(:,2).^2) + pts(:,5).*(pts(:,3).^2) - pts(:,1).*(pts(:,4).^2) +...
        pts(:,4).^3 + pts(:,4).^2 - (1/3)*pts(:,1) + (4/3)*pts(:,4);

return

function vals = magnetism7Vals(pts)
% This method computes evaluations of the 7-variable magnetism polynomial at
% the interpolation points.
% --------------------------------------------------------------------------
% USAGE of "magnetism7Vals"
% vals = magnetism7Vals(pts)
% --------------------------------------------------------------------------
% INPUT
% pts:	interpolation points
%
% OUTPUT
% vals: evaluations of the 7-variable magnetism polynomial at the
%       interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    vals = pts(:,1).^2 + 2*(pts(:,2).^2) + 2*(pts(:,3).^2) + 2*(pts(:,4).^2) +...
        2*(pts(:,5).^2) + 2*(pts(:,6).^2) + 2*(pts(:,7).^2) - pts(:,1);

return

function vals = heartVals(pts)
% This method computes evaluations of the heart dipole polynomial at
% the interpolation points.
% --------------------------------------------------------------------------
% USAGE of "heartVals"
% vals = heartVals(pts)
% --------------------------------------------------------------------------
% INPUT
% pts:	interpolation points
%
% OUTPUT
% vals: evaluations of the heart dipole polynomial at the interpolation points
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    vals = -pts(:,1).*(pts(:,6).^3) + 3*pts(:,1).*pts(:,6).*(pts(:,7).^2) -...
        pts(:,3).*(pts(:,7).^3) + 3*pts(:,3).*pts(:,7).*(pts(:,6).^2) -...
        pts(:,2).*(pts(:,5).^3) + 3*pts(:,2).*pts(:,5).*(pts(:,8).^2) -...
        pts(:,4).*(pts(:,8).^3) + 3*pts(:,4).*pts(:,8).*(pts(:,5).^2) -...
        0.9563453;

return
