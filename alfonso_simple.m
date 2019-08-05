%ALFONSO_SIMPLE is a simple, user-expandable interface for alfonso for
%  solving standard form problems involving pre-defined primitive cones
%
% ADD NEW CONES BY ADDING NEW CASES TO THE FUNCTIONS
%   x0_K, gH_K, AND dimPreprocess (if necessary) BELOW
%
% -------------------------------------------------------------------------
%
% USAGE of "alfonso_simple"
% results = alfonso_simple(c, A, b, K, x0, opts)
% -------------------------------------------------------------------------
% INPUT
% c:                 cost vector
% A:                 constraint matrix
% b:                 right-hand side vector
% K:                 cone, given by a cell array with fields defining the
%                       product of primitive cones to optimize over.
%                    Example:
%                        K{1}.type = 'socp';    % second-order cone
%                        K{1}.dim  = 10;
%                        K{2}.type = 'exp';     % exponential cone, always 3-dimensional
%                        K{3}.type = 'lp';      % nonnegative orthant
%                        K{3}.dim  = 10;
% x0:                initial primal iterate, pass [] to let alfonso choose it
% opts:              algorithmic options, see alfonso.m for details
%
% OUTPUT
% results:           final solution and iteration statistics
%                    see alfonso.m for details
%
% -------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz    <syildiz@email.unc.edu>  
%
% Date: 08/01/2019
% 
% This code has been developed and tested with Matlab R2018a.
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% None.
% --------------------------------------------------------------------------


function results = alfonso_simple(c, A, b, K, x0, opts)

    % preprocess nonnexisting 'dim' fields and 0-dim cones
    [K, Kdims] = dimPreprocess(K);
    
    inputCheck(A, b, c, K);
    
    probData   = struct('c', c, 'A', A, 'b', b, 'Kdims', Kdims);
    
    param.K    = K;
    if isempty(x0)
        x0 = x0_K(K);
    else
        % we should check here what the user gave, but we don't at the moment
    end
    results    = alfonso(probData, x0, @gH_K, param, opts);
    
return


function x0 = x0_K(K)
% x0_K(K) returns a central vector in K: a concatenation of pre-determined
% central points of its component simple cones
%
% adding new cones requires adding a new case below

    n = 0;
    for i = 1:length(K)
        n = n + K{i}.dim;
    end
    
    x0 = zeros(n,1);
    n = 0;
    for i = 1:length(K)
        switch K{i}.type             %%% ADD YOUR CONES HERE (1/3)
            case {'l', 'lp'}
                x = ones(K{i}.dim,1);
            case {'soc', 'socp'}
                x = [1; zeros(K{i}.dim-1,1)];
            case {'rsoc'}
                x = [sqrt(1/2); sqrt(1/2); zeros(K{i}.dim-2,1)];
            otherwise
                error(['unsupported cone type: ', K{i}.type]);
        end
        x0(n+1 : n+K{i}.dim) = x;
        n = n + K{i}.dim;
    end

return

% gH_K(K) implements a barrier function for K computed from component simple cones
%
% INPUT
% x:                 point at which the barrier is to be computed
% params:            optional parameters to use in barrier computation:
%  - params.K        mandatory parameter specifying the cone as a cell array
%
% OUTPUT
% in:                1 (true) if x is in the interior of the cone, 0 (false) otherwise
% g:                 gradient of the barrier function at x;   NaN if in = 0
% H:                 Hessian of the barrier function at x;    NaN if in = 0
% L:                 Cholesky factor of H (lower triangular); NaN if in = 0                    
%
% adding new cones requires adding a new case below
% Tips:
%   - use nargout to only compute whatever is necessary
%   - using chol(H,'lower') for computing L is okay, but one can often do better

function [in, g, H, L] = gH_K(x, params)

    K  = params.K;
    nK = length(K);
    
    in  = true;
    g   = zeros(size(x));
    Hs  = cell(nK,1);  % doesn't do much, but helps with a warning
    Ls  = cell(nK,1);  % doesn't do much, but helps with a warning
    idx = 0; % x subvector index
    for i = 1:nK
        ni = K{i}.dim;
        xi = x(idx+1:idx+ni);
        
        switch K{i}.type              %%% ADD YOUR CONES HERE (2/3)
            case {'l', 'lp'}
                gH      = @gH_LP;
                params  = [];         % the dimension is the only parameter, but that is obtained directly from the length of the vector
            case {'soc', 'socp'}
                gH      = @gH_SOCP;
                params  = [];         % the dimension is the only parameter, but that is obtained directly from the length of the vector
            case {'rsoc'}
                gH      = @gH_RSOC;
                params  = [];         % the dimension is the only parameter, but that is obtained directly from the length of the vector
            %case {'exp'}
            %    gH      = @gH_Exp;   % coming very soon... feel free to implement your own in the meantime :)
            %    params  = [];
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

% Makes sure that every cone in K has a field 'dim' that contains the
% dimension of said cone. Useful for cones where the dimension is either
% not a parameter (e.g., the exponential cone) or is not the natural
% parameter (e.g., the SDP cone, whose dimension is n^2 or n*(n+1)/2,
% depending on the implementation)
%
% may be used to clean up other edge cases, e.g., 0-dimensional cones

function [K, Kdims] = dimPreprocess(K)   %%% ADD YOUR CONES HERE 3/3 (not needed to change this for most cones)

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
    
    N = 0;
    for i = 1:length(K)
        N = N + K{i}.dim;
    end
    if n ~= N
        error('Dimension of (column) vector c must match the dimension of the cone K');
    end
    
return


