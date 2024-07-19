%ALFONSO_SIMPLE is a simple, user-exandable interface for alfonso for
%  solving standard form problems involving pre-defined primitive cones
%
% ADD NEW CONES BY ADDING NEW CASES TO THE FUNCTIONS
%   x0_K, gH_K, and (if necessary) dimPreprocess IN THIS FILE
%
% -------------------------------------------------------------------------
%
% USAGE of "alfonso_simple"
% results = alfonso_simple(c, A, b, K, x0, opts)
% -------------------------------------------------------------------------
% INPUT
% c:                    cost vector
% A:                    constraint matrix
% b:                    right-hand side vector
% K:                    cone, given by a cell array with fields defining
%                          the product of primitive cones to optimize over.
%                       Example:
%                           K{1}.type = 'socp';    % second-order cone
%                           K{1}.dim  = 10;
%                           K{2}.type = 'exp';     % exponential cone, always 3-dimensional
%                           K{3}.type = 'lp';      % nonnegative orthant
%                           K{3}.dim  = 10;
% x0:                   initial primal iterate, pass [] to be chosen by alfonso
% opts:                 algorithmic options, see alfonso.m for details
%   - opts.preprocess:  1 (true) to clean up the input and identify sparse structures,
%                       0 (false) otherwise. Default value: 1
%   - opts.ignorelp:    1 (true) to leave orthants alone in preprocessing.
%                       Default value: 0 (false)
%
% OUTPUT
% results:              final solution and iteration statistics
%                       see alfonso.m for details
%
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


function results = alfonso_simple(c, A, b, K, x0, opts)

    switch nargin
        case 4
            x0 = [];
            opts = struct('preprocess', true);
        case 5
            opts = struct('preprocess', true);
        case 6
            if ~isfield(opts, 'preprocess')
                opts.preprocess = true;
            end
        otherwise
            error('alfonso_simple needs 4-6 input arguments');
    end
    if ~isfield(opts, 'ignorelp')
        opts.ignorelp = false;
    end
            
    % clean up the cone
    if opts.preprocess
        K = dimPreprocess(K);
    end
    
    if isempty(x0)
        x0 = x0_K(K);
    end
    
    inputCheck(A, b, c, K);
    
    % basic preprocessing
    if opts.preprocess
        [c, A, b, K, Kdims, x0, backTrafo] = problemPreprocess(c, A, b, K, x0, opts.ignorelp);
        if issparse(A)
            [rB, frB, As] = sparsePreprocess(A, Kdims);
            probData = struct('c', c, 'A', A, 'b', b, 'Kdims', Kdims, 'rB', rB, 'frB', {frB}, 'As', {As});
            param.sparseA = true;
        else
            probData = struct('c', c, 'A', A, 'b', b, 'Kdims', Kdims);
            param.sparseA = false;
        end
    else
        Kdims = dims(K);
        probData = struct('c', c, 'A', A, 'b', b, 'Kdims', Kdims);
        param.sparseA = false;
    end
    
    % solve the problem with alfonso    
    param.K    = K;
    param.dims = Kdims;
    results    = alfonso(probData, x0, @gH_K, param, opts);
    
    % transform the solution back
    if opts.preprocess
        results.x = results.x(backTrafo);
        results.s = results.s(backTrafo);
    end
    
return


% x0_K(K) returns a central vector in K: a concatenation of pre-determined
% central points of its component simple cones
function x0 = x0_K(K)

    n = 0;
    for i = 1:length(K)
        n = n + K{i}.dim;
    end
    
    x0 = zeros(n,1);
    n = 0;
    for i = 1:length(K)
        switch K{i}.type
            case {'l', 'lp'}
                x = ones(K{i}.dim,1);
            case {'free'}
                x = zeros(K{i}.dim,1);
            case {'soc', 'socp'}
                x = [1; zeros(K{i}.dim-1,1)];
            case {'rsoc'}
                x = [sqrt(1/2); sqrt(1/2); zeros(K{i}.dim-2,1)];
            case {'gpow'}
                x = [ones(length(K{i}.lambda),1); 0];    % not quite sure this is the best choice
            case {'dgpow'}
                x = [ones(length(K{i}.lambda),1); 0];    % not quite sure this is the best choice
            case {'exp'}
                x = [0.7633336255892224; 0.4910129724669193; -0.4197952321239648];  % what else?
            case {'grk1lmi'}
                x = K{i}.x0;
            otherwise
                error(['unsupported cone type: ', K{i}.type]);
        end
        x0(n+1 : n+K{i}.dim) = x;
        n = n + K{i}.dim;
    end

return

% Composite gH function computed from the constituent gH functions. Used as
% the main argument of alfonso().
function [in, g, Hi, Li] = gH_K(x, params)

    K  = params.K;
    nK = length(K);
    
    in  = true;
    g   = zeros(size(x));
    Li  = cell(nK,1);  % doesn't do much, but helps with a warning
    Hi  = cell(nK,1);  % doesn't do much, but helps with a warning
    
    idx  = 0;           % x subvector index
    for i = 1:nK
        ni = K{i}.dim;
        xi = x(idx+1:idx+ni);

        switch K{i}.type
            case {'l', 'lp'}
                gH    = @gH_LP;
                par0  = [];    % the dimension is the only parameter, but that is obtained directly from the length of the vector
            case {'soc', 'socp'}
                gH    = @gH_SOCP;
                par0  = [];    % the dimension is the only parameter, but that is obtained directly from the length of the vector
            case {'rsoc'}
                gH    = @gH_RSOC;
                par0  = [];    % the dimension is the only parameter, but that is obtained directly from the length of the vector
            case {'gpow'}
                gH    = @gH_GPow;
                par0  = K{i}.lambda;
            case {'dgpow'}
                gH    = @gH_DGPow;
                par0  = K{i}.lambda;
            case {'exp'}
                gH    = @gH_Exp;
                par0  = [];
            case {'grk1lmi'}
                    gH    = @gH_rk1LMI;
                    par0  = K{i}; % struct('Vs', K{i}.Vs, 'ws', K{i}.ws, 'nonneg', K{i}.nonneg, 'ext', K{i}.ext);
            otherwise
                error(['unsupported cone type: ', K{i}.type]);
        end

        switch nargout
            case 4
                [in0,g0,Hi0,Li0] = gH(xi, par0);
            case 3
                [in0,g0,Hi0] = gH(xi, par0);
            case {1,2}
                [in0,g0] = gH(xi, par0);
        end

        if in0
            g(idx+1:idx+ni) = g0;
            if nargout >= 3
                Hi{i} = Hi0;
                if nargout == 4
                    Li{i} = Li0;
                end
            end
        else
            in  = 0;
            g   = NaN;
            Hi  = NaN;
            Li  = NaN;
            return;
        end

        idx = idx+ni;
    end

    if ~params.sparseA
        if nargout > 2
            Hi = @(v)(concatF(Hi,params.dims,v));
            if nargout > 3
                Li = @(M)(concatF(Li,params.dims,M));
            end
        end
    end
    
return

function FofM = concatF(Fs, Kdims, M)

    if issparse(M)
        
        FMs = cell(length(Kdims),1);
        idx = 0;  % x subvector index 
        for i=1:length(Kdims)
            FMs{i} = Fs{i}(M(idx+1:idx+Kdims(i), :));
            idx = idx + Kdims(i);
        end
        FofM = vertcat(FMs{:});
        
    else
        
        FofM = zeros(size(M));
        idx = 0;  % x subvector index
        for i=1:length(Kdims)
            FofM(idx+1:idx+Kdims(i), :) = Fs{i}(M(idx+1:idx+Kdims(i), :));
            idx = idx + Kdims(i);
        end
        
    end
    
return


function K = dimPreprocess(K)
% This method performs some basic preprocessing on the cone structure K for
% efficiency and compatibility with other functions. It should NOT be
% called outside of alfonso_simple().
%
% Included functionality:
%  - add missing 'dim' fields to cones for which the dimension is not the
%    natural parameter (e.g., exponential cone, semidefinite cone)
%  - remove 0-dimensional cones
%
% --------------------------------------------------------------------------
% USAGE of "dimPreprocess"
% K = dimPreprocess(K)
% --------------------------------------------------------------------------
% INPUT
% K:            cone structure
%
% OUTPUT
% K:            updated cone structure
% --------------------------------------------------------------------------

    for i=length(K):-1:1
        if ~isfield(K{i},'dim')
            switch K{i}.type
                case 'exp'
                    K{i}.dim = 3;
                case 'grk1lmi'
                    if K{i}.ext
                        K{i}.dim = size(K{i}.ws,1) + 1;
                    else
                        K{i}.dim = size(K{i}.ws,1);
                    end
                otherwise
                    error(['missing dimension in cone ', int2str(i), ', type: ', K{i}.type]);
            end
        end
        if K{i}.dim <= 0
            K(i) = [];
        end
    end
   
return

function [c, A, b, K, Kdims, x0, backTrafo] = problemPreprocess(c, A, b, K, x0, ignoreLP)
% This method performs some basic preprocessing on the primal problem for
% efficiency and compatibility with other functions. It should NOT be
% called outside of alfonso_simple().
%
% Included functionality:
%  - move free variables into a single second-order cone
%  - move nonnegative variables into a single orthant
%
% --------------------------------------------------------------------------
% USAGE of "problemPreprocess"
% [c, A, b, K, Kdims, x0, backTrafo] = problemPreprocess(c, A, b, K, x0, ignore)
% --------------------------------------------------------------------------
% INPUT
% c, A, b, K:   problem data and cone structure
% x0:           initial point
% ignoreLP:     if true (1), it does not mess with the orthants
%
% OUTPUT
% c, A, b, K:   updated problem data and cone structure
% x0:           updated initial point
% Kdims:        array of integers, contains the dimensions of K{i}
% backTrafo:    variable permutation to obtain the solution of the original
%               problem from the solution of the preprocessed problem
% --------------------------------------------------------------------------

    [m,n] = size(A);
    freeVars   = false(n,1);
    nonnegVars = false(n,1);
    keepCones  = false(length(K),1);
    
    % find all the variables we might want to move
    idx = 0;
    for i=1:length(K)
        if strcmpi(K{i}.type, 'free')
            freeVars(idx+1 : idx+K{i}.dim) = true;
        elseif ~ignoreLP && (strcmpi(K{i}.type, 'l') || strcmpi(K{i}.type, 'lp'))
            nonnegVars(idx+1 : idx+K{i}.dim) = true;
        else
            keepCones(i) = true;
        end
        idx = idx + K{i}.dim;
    end
    
    otherVars = ~(freeVars | nonnegVars);
    dummy     = any(freeVars);
    
    % new variable order: dummy variable, free variables, nonnegatives, rest
    c  = [zeros(1,dummy);  c(freeVars);    c(nonnegVars);    c(otherVars)];
    A  = [zeros(m,dummy),  A(:,freeVars),  A(:,nonnegVars),  A(:,otherVars)];
    
    xf = x0(freeVars) ;
    x0 = [zeros(1,dummy); xf; x0(nonnegVars);   x0(otherVars)];
    if dummy
        x0(1) = 1 + norm(xf);
    end
    
    K0 = {};
    if any(freeVars)
        K0 = horzcat(K0, struct('type', 'socp', 'dim', sum(freeVars)+1));
    end
    if any(nonnegVars)
        K0 = horzcat(K0, struct('type', 'lp', 'dim', sum(nonnegVars)));
    end
    K = horzcat(K0, K(keepCones));
    
    % compute the permutation of variables that arranges everything back
    ind = 1:n;
    ind = [ind(freeVars), ind(nonnegVars), ind(otherVars)];
    if dummy
        ind = [n+1, ind];
        ind(ind) = 1:n+1;
        ind = ind(1:end-1);
    else
        ind(ind) = 1:n;
    end
    backTrafo = ind;
    
    Kdims = dims(K);
    
return

% Identifies the block structure of A.
% The As output is optional in case we don't really want to store the dense
% components??
% As: dense components of A' (transpose!)
% rB: "row blocks", logical array (number of rows of A x number of blocks)
% frB: columnwise find(rB)

function [rB, frB, As] = sparsePreprocess(A, Kdims)

    nK = length(Kdims);
    if nargout == 3
        As = cell(nK,1);
    end

    rB  = zeros(size(A,1), nK, 'logical');
    frB = cell(nK,1);
    sidx = 1;                       % starting index of the column block
    for i=1:nK
        eidx = sidx - 1 + Kdims(i); % ending index of the column block
        rB(:,i) = any(A(:,sidx:eidx), 2);
        frB{i} = find(rB(:,i));
        if nargout == 3
            As{i} = full(A(frB{i},sidx:eidx).');
        end
        sidx = eidx + 1;
    end

return


function Kdims = dims(K)

    Kdims = zeros(length(K),1);
    for i=1:length(K)
        Kdims(i) = K{i}.dim;
    end
    
return


function inputCheck(A, b, c, K)

    [m, n] = size(A);

    if m <= 0 || n <= 0
        error('Input matrix A must be nontrivial.');
    end
    if m ~= size(b,1)
        error('Dimension of (column) vector b must match the number of rows in A.');
    end
    if n ~= size(c,1)
        error('Dimension of (column) vector c must match the number of columns in A.');
    end

    if n ~= sum(dims(K))
        error('Dimension of (column) vector c must match the dimension of the cone K');
    end
    
return


