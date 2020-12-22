% This code formulates and solves a portfolio optimization with a
% factor risk model and market impact model, using the simple interface of
% alfonso. See the description of the model in the header of PORTFOLIO(f,n).
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% None.
% -------------------------------------------------------------------------

function results = portfolio(f, n)
% This is the main method for the portfolio optimization problem.
% -------------------------------------------------------------------------
% USAGE of "portfolio"
% results = portfolio(f, n)
% -------------------------------------------------------------------------
% INPUT
% f:                    number of factors
% n:                    number of assets
%
% OUTPUT
% results:              final solution and iteration statistics
%                       see alfonso.m for details
%
% BRIEF MODEL DESCRIPTION:
%
% The asset return covariance matrix has the factor model structure
%   Sigma = B Omega B' + Delta,
% where B is an n x f matrix, Omega is an f x f PSD matrix, and Delta is a
% diagonal positive definite  matrix.
%
% Factor risk model parameters:
%  - B:         factor exposure matrix
%  - Omega:     factor covariance matrix
%  - diagDelta: specific variance of asset returns
% 
% Market impact model parameters:
%  - beta:      exponent
%  - lambda:    multipliers
%
% Additional model parameters:
%  - alpha: asset expected returns (alpha's)
%  - portfolio risk limit: gamma
%  - initial portfolio: h
%
% Decision variables:
%  - optimal portfolio weights: x
%  - transaction weights: t
%
% The optimization model is of the form
%
%   maximize_{x,t}   alpha'*x - sum_i (lambda_i*|t|^beta)
%   subject to       sum(x) == 1
%                    x >= 0
%                    x+t = h
%                    sqrt(x'*(B*Omega*B' + Delta*x) <= gamma
%                    
% -------------------------------------------------------------------------


    % generate some realistic-looking random data
    [alpha, B, Omega, diagDelta, beta, lambda, gamma, h] = getRandomData(f,n);
    
    %[n, f] = size(B);
    
    % variable order: x, (z_i, 1, t_i)_{i=1:n}, (gamma, d, u)
    
    % K = R_+^n x P_{1/beta,1-1/beta}^n x Q_{n+f+1}
    K = cell(1,n+2);
    K{1} = struct('type','lp','dim',n);
    [K{2:n+1}] = deal(struct('type','gpow','dim',3,'lambda',[1/beta; 1-1/beta]));
    K{n+2} = struct('type','socp','dim',n+f+1);
    
    rows = 1+n+n+f+n+1;
    cols = n+3*n+(n+f+1);
    
    A = zeros(rows, cols);
    b = zeros(rows,1);
    c = zeros(cols,1);
    
    c(1:n) = -alpha;
    c(n+1:3:4*n-2) = lambda;
    
    % e'x = 1
    A(1,1:n) = 1;
    b(1) = 1;

    % x-t = h
    A(2:n+1,1:n) = eye(n);
    A(2:n+1,n+3:3:n+3*n) = -eye(n);
    b(2:n+1) = h;

    % d - Delta^{1/2}*x = 0
    A(n+2:2*n+1,4*n+1+(1:n)) = eye(n);
    A(n+2:2*n+1,1:n) = -diag(sqrt(diagDelta));

    % u - Omega^{1/2}B'*x = 0
    A(2*n+2:2*n+1+f,5*n+1+(1:f)) = eye(f);
    A(2*n+2:2*n+1+f,1:n) = -sqrtm(Omega)*B';

    % constants
    A(2*n+1+f+(1:n), n+2:3:4*n-1) = eye(n);
    b(2*n+1+f+(1:n)) = 1;
    A(2*n+1+f+n+1,4*n+1) = 1;
    b(2*n+1+f+n+1,1) = gamma;
    
    fprintf('A size: %dx%d;  sparsity:%f\n',rows,cols,nnz(A)/numel(A));
    
    A = sparse(A);
    
    results = alfonso_simple(c, A, b, K);
    results.data = struct('h',h,'alpha',alpha,'beta',beta,'gamma',gamma,'Delta',sparse(diag(diagDelta)),'B',B,'Omega',Omega);

return

function [alpha, B, Omega, diagDelta, beta, lambda, gamma, h] = getRandomData(f,n)
% Utility method to generate random data for the portfolio optimization problem.
% -------------------------------------------------------------------------
% USAGE of "getRandomData"
% [alpha, B, Omega, diagDelta, beta, lambda, gamma, h] = getRandomData(f, n)
% -------------------------------------------------------------------------

    seed = 2020;
    if isOctave()
        randn("state", seed);
    else
        rng(seed,'twister');        
    end

    % market impact model parameters
    lambda = 0.1*ones(n,1);        % multipliers
    beta  = 5/3;                   % exponent
    
    % asset expected returns (alpha's)
    alpha = randn(n,1);
    
    % risk model parameters
    % asset return covariance matrix has the factor model structure
    %   B Omega B' + Delta,
    % where B is an n x f matrix, Omega is an f x f PSD matrix, and Delta
    % is a diagonal PSD matrix
    
    B  = [ones(n,1), 0.5*randn(1,f-1) + 0.2*randn(n,f-1)];  % first column represents the market intercept
    
    % generate a random orthogonal matrix Q from the Haar distribution
    [Q,R] = qr(randn(f));
    Q = Q*diag(sign(diag(R)));
    
    Omega = Q * diag((0.15 + 0.03*randn(f,1)).^2) * Q';
    diagDelta = (0.2 + 0.05*randn(n,1)).^2; % Delta = diag(diagDelta)
    
    % optimization model parameters
    h     = 1/n * ones(n,1);   % equally weighted initial portfolio
    Sigma = B*Omega*B' + diag(diagDelta);
    gamma = sqrt(h'*Sigma*h);  % portfolio risk limit = current risk

    disp(['Portfolio risk limit (gamma) = ', num2str(gamma)]);
return
