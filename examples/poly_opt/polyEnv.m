% This code is an implementation of the sum-of-squares optimization approach 
% based on non-symmetric conic optimization and polynomial interpolants 
% presented in:
%
% D. Papp and S. Yildiz. Sum-of-squares optimization without semidefinite 
% programming. Available at https://arxiv.org/abs/1712.01792.
%
% The implementation formulates and solves the polynomial envelope problem
% described in the same reference.
%
% Note: the implementation follows the paper above. With the current
% version of alfonso, using alfonso_simple and the `rk1LMI' cone is
% recommended, as it is far simpler and more efficient.
% -------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz
%
% Date: 2024/07/15
%
% This code has been developed and tested with Matlab R2023b.
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% None.
% -------------------------------------------------------------------------


function results = polyEnv(intParams, numPolys, degPolys, tol, seed)
% This is the main method for the sum-of-squares optimization approach to 
% to the polynomial envelope problem.
% --------------------------------------------------------------------------
% USAGE of "polyEnv"
% results = polyEnv(intParams, numPolys, degPolys, tol, seed)
% Use seed = 2017 to reproduce the results in the reference above.
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
% numPolys:         number of approximated polynomials
% degPolys:         degree of approximated polynomials. it is assumed that
%                   degPolys <= intParams.d.
% tol:              tolerance parameter for optimization
% seed:             seed for the random number generator
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

    n   = intParams.n;
    d   = intParams.d;
    U   = intParams.U;
    L   = intParams.L;
    P   = intParams.P;
    pts = intParams.pts;
    
    % ORDER OF VARIABLES 
    % x = [x_1; x_2; ...] \in WSOS_(n,2*d)^* x WSOS_(n,2*d)^* x ...
    % x_1 corresponds to the 1st approximated polynomial, x_2 to the 2nd,...
    
    % WEIGHTS
    % weight #0: 1
    % degree of associated SOS multiplier: 2*d
    % weight #j for j = 1,...,n: (1-t_j^2)
    % degree of associated SOS multiplier: 2*d-2
    
    LWts = repmat(nchoosek(n+d-1,n),n,1);
    
    % PARAMETERS ASSOCIATED WITH THE METHOD gH_polyEnv
    % see input description for gH_polyEnv for details
    gH_Params.n = n;
    gH_Params.d = d;
    gH_Params.U = U;
    gH_Params.numPolys = numPolys;
    
    gH_Params.L     = L;
    gH_Params.LWts  = LWts;
    nu              = numPolys*(L+sum(LWts));
    gH_Params.bnu	= nu+1;
    
    % evaluations of the weights 1-t_j^2 at the interpolation points:
    % the first column is for the weight 1-t_1^2, the second is for
    % 1-t_2^2,...
    wtVals = 1-pts.^2;

    % P has orthonormal columns:
    gH_Params.P = P;
    % associated positive semidefinite cone constraints:
    % P'*diag(x_1)*P >= 0,
    % P'*diag(x_2)*P >= 0,...
    PWts = cell(n,1);
    for j = 1:n
        PWts{j}         = diag(sqrt(wtVals(:,j)))*P(:,1:LWts(j));
        [PWts{j}, ~]    = qr(PWts{j}, 0);
        % associated positive semidefinite cone constraints: 
        % PWts{j}'*diag(x_1)*PWts{j} >= 0,
        % PWts{j}'*diag(x_2)*PWts{j} >= 0,...
    end
    gH_Params.PWts = PWts;
    
    % x_1 \in WSOS_(n,2*d)^* <=> 
    % P'*diag(x_1)*P >= 0, PWts{1}'*diag(x_1)*PWts{1} >= 0, PWts{2}'*diag(x_1)*PWts{2} >= 0,...
    % x_2 \in WSOS_(n,2*d)^* <=> 
    % P'*diag(x_1)*P >= 0, PWts{1}'*diag(x_2)*PWts{1} >= 0, PWts{2}'*diag(x_2)*PWts{2} >= 0,...
    
    % DATA FOR THE CONIC OPTIMIZATION PROBLEM
    probData.A = sparse(repmat(eye(U),1,numPolys));
    probData.b = intParams.w;
    probData.c = genRandPolyVals(n, numPolys, degPolys, intParams.P0);
      
    tic
    
    % INITIAL PRIMAL ITERATE
    x0 = ones(numPolys*U,1);
    [~, g0, ~, ~] = gH_polyEnv(x0, gH_Params);
    % scaling factor for the primal problem
    rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
    % scaling factor for the dual problem
    rD = max((1+abs(g0))./(1+abs(probData.c)));
    % initial primal iterate
    x0 = repmat(sqrt(rP*rD),numPolys*U,1);  
    
    % CUSTOM ALGORITHMIC OPTIONS
    opts.optimTol = tol;
    
    % CALL TO alfonso 
    results = alfonso(probData, x0, @gH_polyEnv, gH_Params, opts);
    fprintf('alfonso is done.\n');
    
    toc
    
    % prints error statistics
    fprintf('\n');
    fprintf('FINAL:\n');
    fprintf('Relative primal infeasibility: %d\n', results.rel_pIn);
    fprintf('Relative dual infeasibility: %d\n', results.rel_dIn);
    fprintf('Relative duality gap: %d\n', results.rel_dGap);
    fprintf('Relative complementarity gap: %d\n\n', results.rel_cGap);
    
return

function [in, g, Hi, Li] = gH_polyEnv(x, params)
% This method computes the gradient and Hessian of the barrier function for
% the problem of polynomial envelopes.
% --------------------------------------------------------------------------
% USAGE of "gH_polyEnv"
% [in, g, H, L] = gH_polyEnv(x, params)
% --------------------------------------------------------------------------
% INPUT
% x:                    primal iterate
% params:               parameters associated with the method gH
% - params.n:           number of arguments to the polynomials
% - params.d:           largest degree of polynomials to be squared
% - params.U:           dimension of the space of (params.n)-variate 
%                       degree-(2*params.d) polynomials
% - params.numPolys:    number of approximated polynomials
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
%                       params.PWts{j} is sqrt(1-t_j^2).
%
% OUTPUT
% in:	0 if x is not in the interior of the cone. 1 if x is in the 
%       interior of the cone.
% g:    gradient of the barrier function at x
% Hi:	function representing the inverse Hessian action at x
% Li:   function representing the inverse Cholesky action or similar
% (see the built-in barrier functions or the documentation for more info)
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    in = 1;
    n = params.n;
    U = params.U;
    numPolys = params.numPolys;
    P = params.P;
    PWts = params.PWts;
    
    g = zeros(numPolys*U,1);
    L = cell(numPolys,1);

    % ORDER OF VARIABLES 
    % x = [x_1; x_2; ...] \in WSOS_(n,2*d)^* x WSOS_(n,2*d)^* x ...
    % x_1 corresponds to the 1st approximated polynomial, x_2 to the 2nd,...
    
    off = 0;
    for polyId = 1:numPolys
        xPoly = x(off+(1:U));
        % for the weight 1
        [inPoly, gPoly, HPoly] = gH_SOSWt(xPoly,P);

        if inPoly == 1
            for j = 1:n
                % for the weight 1-t_j^2
                [inPolyWt, gPolyWt, HPolyWt] = gH_SOSWt(xPoly,PWts{j});
                inPoly  = inPoly & inPolyWt;
                if inPoly == 1
                    gPoly   = gPoly+gPolyWt;
                    HPoly   = HPoly+HPolyWt;
                else
                    gPoly   = NaN;
                    HPoly   = NaN;
                    break;
                end
            end
        end

        if inPoly == 1
            % checks positive semidefiniteness of HPoly one more time.
            % HPoly may fail Cholesky factorization even if inPoly == 1
            % due to numerical errors in summing up HPolyWt's above.
            [LPoly, err] = chol(HPoly,'lower');
            inPoly  = inPoly & (err == 0);
        end

        if inPoly == 1
            g(off+(1:U)) = gPoly;
            L{polyId} = LPoly;
            off = off + U;
        else
            in  = 0;
            g  = NaN;
            Li = NaN;
            Hi = NaN;
            return;
        end
    end

    if nargout >= 3
        Hi = @(v)(concatH(L,v));
        if nargout == 4
            Li = @(M)(concatL(L,M));
        end
    end

return

function LiM = concatL(Ls, M)
    
    U = size(Ls{1},1);
    
    LiMs = cell(length(Ls),1);
    idx = 0;  % x subvector index
    for i=1:length(Ls)
        LiMs{i} = Ls{i} \ M(idx+1:idx+U, :);
        idx = idx + U;
    end
    LiM = vertcat(LiMs{:});

return

function Hiv = concatH(Ls, v)

    U = size(Ls{1},1);
    
    Hiv = zeros(size(v));
    idx = 0;  % x subvector index
    for i=1:length(Ls)
        Hiv(idx+1:idx+U) = Ls{i}'\(Ls{i}\v(idx+1:idx+U));
        idx = idx + U;
    end

return



function [in, g, H] = gH_SOSWt(x, P)
% This method computes the gradient and Hessian of a barrier function term
% corresponding to a single weight and single approximated polynomial for
% the problem of polynomial envelopes.
% --------------------------------------------------------------------------
% USAGE of "gH_SOSWt"
% [in, g, H] = gH_SOSWt(x, P)
% --------------------------------------------------------------------------
% INPUT
% x:    subvector of the primal iterate corresponding to a single
%       approximated polynomial
% P:    evaluations of "weighted" basis polynomials at the
%       interpolation points
%
% OUTPUT
% in:   0 if P'*diag(x)*P is not positive definite. 1 if P'*diag(x)*P is
%       positive definite.
% g:    gradient of the barrier function term corresponding to a single
%       weight and single approximated polynomial at x
% H:    Hessian of the barrier function term corresponding to a single
%       weight and single approximated polynomial at x
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    Y = P'*diag(x)*P;

    if ~issymmetric(Y)
        Y = (Y+Y')/2;
    end
    
    [L, err] = chol(Y, 'lower');
    if err > 0
        in = 0;
        g = NaN;
        H = NaN;
    else
        in = 1;
        V = L\P';
        H = V'*V;  % VtV; not the Hessian yet, but below we compute the Hessian in-place from this.

        g = -diag(H);
        H = H.^2;
    end
    
return

function vals = genRandPolyVals(n, numPolys, degPolys, P)
% This method computes evaluations of randomly generated polynomials at
% the interpolation points.
% --------------------------------------------------------------------------
% USAGE of "genRandPolyVals"
% vals = genRandPolyVals(n, numPolys, degPolys, P)
% --------------------------------------------------------------------------
% INPUT
% n:            number of arguments to the polynomials
% numPolys:     number of polynomials to be generated
% degPolys:     degree of polynomials to be generated. it is assumed that
%               degPolys <= d.
% P:            evaluations of a basis for the space of n-variate degree-d 
%               polynomials at the interpolation points. it is assumed that 
%               the first nchoosek(n+degPolys,degPolys) columns are a basis
%               for the space of n-variate degree-(degPolys) polynomials.
%
% OUTPUT
% vals: evaluations of numPolys randomly generated n-variate
%       degree-(degPolys) polynomials at the interpolation points.
%       the evaluation vectors are stacked vertically. 
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    % generates coefficients in the basis that indexes the columns of P
    LDegs = nchoosek(n+degPolys, n);
    coeffs = randi([-9,9], LDegs, numPolys);
    vals = P(:,1:LDegs)*coeffs;
    
    vals = vals(:);
    
return
