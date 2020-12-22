%ALFONSO_SIMPLE is a simple, user-extensible interface for alfonso for
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
%   - opts.preprocess:  1 (true) to clean up the input, 0 (false) otherwise.
%                       Default value: 1
%   All other options are passed to alfonso.
% -------------------------------------------------------------------------
% OUTPUT
% results:              final solution and iteration statistics
%                       see alfonso.m for details
%
% -------------------------------------------------------------------------
%
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
% Version: 07/20/2020
% 
% This code has been developed and tested with Matlab R2018a.
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% None.
% --------------------------------------------------------------------------


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
        [c, A, b, K, Kdims, x0, backTrafo] = problemPreprocess(c, A, b, K, x0);
    else
        Kdims = dims(K);
    end
    
    % solve the problem with alfonso
    probData   = struct('c', c, 'A', A, 'b', b, 'Kdims', Kdims);
    param.K    = K;
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
                x = [ones(length(K{i}.lambda),1); 0]; % David's favorite
                %x = [sqrt(1+K{i}.lambda); 0];        % more standard
            case {'dgpow'}
                x = [ones(length(K{i}.lambda),1); 0]; % David's favorite
                %x = [sqrt(1+K{i}.lambda); 0];        % more standard
            case {'exp'}
                x = [0.7633336255892224; 0.4910129724669193; -0.4197952321239648];  % David's favorite
                %x = [1.290927709856958; 0.8051020015847954; -0.8278383990656786];  % more standard
            otherwise
                error(['unsupported cone type: ', K{i}.type]);
        end
        x0(n+1 : n+K{i}.dim) = x;
        n = n + K{i}.dim;
    end

return


function [in, g, H, L] = gH_K(x, params)

    K  = params.K;
    nK = length(K);
    
    in  = true;
    g   = zeros(size(x));
    Hs  = cell(nK,1);  % doesn't do much, but helps with a warning
    Ls  = cell(nK,1);  % doesn't do much, but helps with a warning
    idx = 0;           % x subvector index
    for i = 1:nK
        ni = K{i}.dim;
        xi = x(idx+1:idx+ni);
        
        switch K{i}.type
            case {'l', 'lp'}
                gH      = @gH_LP;
                params  = [];         % the dimension is the only parameter, but that is obtained directly from the length of the vector
            case {'soc', 'socp'}
                gH      = @gH_SOCP;
                params  = [];         % the dimension is the only parameter, but that is obtained directly from the length of the vector
            case {'rsoc'}
                gH      = @gH_RSOC;
                params  = [];         % the dimension is the only parameter, but that is obtained directly from the length of the vector
            case {'gpow'}
                gH      = @gH_GPow;
                params  = K{i}.lambda;
            case {'dgpow'}
                gH      = @gH_DGPow;
                params  = K{i}.lambda;
            case {'exp'}
                gH      = @gH_Exp;
                params  = [];         % no parameters
            otherwise
                error(['unsupported cone type: ', K{i}.type]);
        end
        
        switch nargout
            case 4
                [in0,g0,H0,L0] = gH(xi, params);
            case 3
                [in0,g0,L0]    = gH(xi, params);
            case {1,2}
                [in0,g0]       = gH(xi, params);
        end
        
        if in0
            g(idx+1:idx+ni) = g0;
            if nargout > 2
                Hs{i} = sparse(H0);
                if nargout > 3
                    Ls{i} = sparse(L0);
                end
            end
        else
            in = 0;
            g  = NaN;
            H  = NaN;
            L  = NaN;
            return;
        end
        
        idx = idx+ni;
    end
    
    if nargout > 2
        H = blkdiag(Hs{:});
    end
    
    if nargout > 3
        L = blkdiag(Ls{:});
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
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    for i=length(K):-1:1
        if ~isfield(K{i},'dim')
            switch K{i}.type
                case 'exp'
                    K{i}.dim = 3;
                otherwise
                    error(['missing dimension in cone ', int2str(i), ', type: ', K{i}.type]);
            end
        end
        if K{i}.dim <= 0
            K(i) = [];
        end
    end
    
return

function [c, A, b, K, Kdims, x0, backTrafo] = problemPreprocess(c, A, b, K, x0)
% This method performs some basic preprocessing on the primal problem for
% efficiency and compatibility with other functions. It should NOT be
% called outside of alfonso_simple().
%
% Included functionality:
%  - move free variables into a single second-order cone
%  - move nonnegative variables into a single orthant
%
% --------------------------------------------------------------------------
% USAGE of "dimPreprocess"
% [c, A, b, K, Kdims, x0, backTrafo] = problemPreprocess(c, A, b, K, x0)
% --------------------------------------------------------------------------
% INPUT
% c, A, b, K:   problem data and cone structure
% x0:           initial point
%
% OUTPUT
% c, A, b, K:   updated problem data and cone structure
% x0:           updated initial point
% Kdims:        array of integers, contains the dimensions of K{i}
% trafo
% backTrafo:      
%
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
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
        elseif strcmpi(K{i}.type, 'l') || strcmpi(K{i}.type, 'lp')
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


