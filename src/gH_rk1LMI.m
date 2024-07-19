%GH_RK1LMI implements a membership and barrier function oracle for linear
% matrix inequalities (LMIs), in the variables x or (z,x), of the form
%     V_i*diag(w_i \circ x)*V_i' \succcurlyeq 0 \forall i = 1, ..., k,
% or
%     V_i*diag(w_i \circ x)*V_i' \succcurlyeq z*I_m \forall i = 1, ..., k,
% where each V_i is an m \times n_i matrix and each w_i is an n-vector.
% Equivalently, the latter can be written as
%     \sum_{j=1}^n x(j) * w_i(j) * V(:,j)*V(:,j)^T \succcurlyeq z*I_m
%
% Optionally, x >= 0 may also be required.
%
% The extended cone and the sign constraints are toggled using the input
% parameters (see INPUT below).
%
% Notes:
%   1) Order of variables in the extended cone is [z,x] \in R x R^n.
%   2) If each V_i is full row-rank and w>0, then x = ones(n,1) is strictly
%      feasible with z=0. In general, z=0 is not feasible, but for every x
%      there is a z that works.
%   3) Assumes each function x -> V_i*diag(w \circ x)*V_i' is injective.
%
% --------------------------------------------------------------------------
% USAGE of "gH_rk1LMI"
% [in, g, Hi, Li] = gH_rk1LMI(x, param)
% --------------------------------------------------------------------------
% INPUT
%     x:             primal iterate
%     param.ws:      a matrix whose ith column is the vector w_i above
%     param.Vs:      cell array with the matrices V_i above
%     param.nonneg:  true (1) if x >= 0 is also required, false (0) if x is
%                    not sign-constrained
%     param.ext:     true (1) if the right-hand side has z*I, false (0) if
%                    there is no z variable, and the right-hand side is 0
%
% OUTPUT
% in:       0 if x is not in the interior of the cone. 1 if x is in the
%           interior of the cone.
% g:        gradient of the barrier function at x
% Hi:	    function representing the inverse Hessian action at x
% Li:       function representing the inverse Cholesky action or similar
%
% The last 3 output may be anything if in==0.
% --------------------------------------------------------------------------
% Details on the last two outputs: denoting the Hessian at x by H,
%  - Li(M) returns L\M, where L is any matrix satisfying LL' = inv(H)
% Needs to work for a matrix.
%  - Hi(v) returns H\v.
% Here, v can be assumed to be a column vector (not a matrix).
%
% It is often sensible, but not always necessary or most efficient, to
% implement Hi using Li or an explicitly computed Cholesky factor.
% --------------------------------------------------------------------------
% Copyright (C) 2024 David Papp  <dpapp@ncsu.edu>
% --------------------------------------------------------------------------

function [in, g, Hi, Li] = gH_rk1LMI(zx, param)

% The first input is [z;x] or x.
if param.ext
    z = zx(1);
    x = zx(2:end);
else
    x = zx;
end

g = NaN;  Hi = NaN;  Li = NaN; % will be overwritten if in==1

% In the nonnegative cone?
if param.nonneg && any(x <= 0)
    in = 0; return;
end

[n, numVs] = size(param.ws);

% Initialize g and H (if needed).
% For simplicity, we treat the dimension as (n+1), and will drop the first
% component later if needed.
if nargout > 1
    if param.nonneg
        g = [0; -1./x];
    else
        g = zeros(n+1,1);         % +1 for the z component
    end
    if nargout > 2
        % Variables to accumulate each block of the Hessian.
        H11   = 0;
        H1end = zeros(n,1);
        if param.nonneg
            % Variable to accumulate the Hessian.
            H2end = diag(x.^(-2));
        else
            H2end = zeros(n,n);
        end
    end
end

% THE REST IS ALL ABOUT THE LMIs.

in = 1; % innocent until proven guilty

% barrier = sum of individual barriers; accumulate g and H in for i=... loop
for i=1:numVs
    % compute M := V_i*diag(w_i \circ x)*V_i' - z*eye(n)
    V = param.Vs{i};
    w = param.ws(:,i);
    
    [m,n] = size(V);
    M = V*diag(w.*x)*V.'; % not the final M yet!

    % compute M without vectorization
    % Vx = zeros(m,n);
    % for u = 1:n
    %     Vx(:,u) = w(u)*x(u)*V(:,u);
    % end
    % M = Vx*V'; % not the final M yet! Also, not symmetric!
    % Vx = zeros(m,n);
    % for u = 1:n
    %     Vx(:,u) = sqrt(w(u)*x(u))*V(:,u);
    % end
    % M = Vx*Vx.'; % not the final M yet!

    % -zI?
    if param.ext
        for u=1:m
            M(u,u) = M(u,u) - z;
        end
    end
    % M is now computed.
    
    % factor M to test feasibility
    [L,err] = chol(M,'lower');
    if err
        in = 0;
        return;
    end
    
    if nargout == 1 && (i==numVs)
        % if we are here, we still have in = 1
        return;
    end
    
    if nargout > 1
        
        LinvV = L \ V;
        if param.ext
            Linv   = inv(L);
            g(1) = g(1) + Linv(:)'*Linv(:);  % this is surprising; any way to avoid explicit inversion?
                                             % tr(M^{-1}) = tr(L^{-T}*L^{-1}) = sum of squared entries of L^{-1}
        end
        if nargout == 2
            % accumulate gradient without the remaining Hessian elements
            g(2:end) = g(2:end) - w .* vecnorm(LinvV,2,1)'.^2;
            if (i==numVs)
                if ~param.ext
                    g = g(2:end);
                end
                return;
            end
        end
        
        if nargout >= 3
            
            % First (z) column of the Hessian.
            if param.ext
                % H_zz: multivariate Faà di Bruno: \nabla_{zz} f(z,x) =
                %                                       = tr(M^{-2}) = sum of squared entries of M^{-1}
                invM = Linv.'*Linv;
                H11 = H11 + invM(:)'*invM(:);
                
                % H_{zx}: even uglier calculus vi' * M^{-2} * vi
                MinvV = L' \ LinvV;
                for u=1:n
                    H1end(u) = H1end(u,1) - w(u) * MinvV(:,u)' * MinvV(:,u);
                end
            end
            
            % H_{xx}:
            % gradient wrt x is computed from Hessian diagonal
            Hi             = (LinvV.'*LinvV);                   % not the Hessian yet! only V.'*(M\V) but we need the diag
            g(2:end)       = g(2:end) - w .* diag(Hi);          % extract diag for gradient

            % Hessian wrt x.
            H2end = H2end + (w*w.').*(Hi.^2);
            %H2end = H2end + w.' .* Hi.^2 .* w; % a bit faster, not automatically symmetric.

            if (i==numVs)
                if param.ext
                    H = [H11, H1end'; H1end, H2end];
                else
                    g = g(2:end);
                    H = H2end;
                end
                
                %Hi = @(v)(H\v);

                %Hi = decomposition(H);
                %Hi = @(v) Hi\v;
                
                %[L,D,P] = ldl(H);
                %Hi = @(v)( P*(L'\(D\(L\(P'*v)))) );

                % direct solve with Cholesky
                [L,err] = chol(H,'lower');
                if err
                    in = 0;
                    return;
                else
                    Hi = @(v)(L'\(L\v));
                    if nargout == 4
                       Li = @(M) L\M;
                    end % if nargout == 4
                    return;
                end
            end
            
        end % if nargout >= 3
    end % if nargout > 1
end % for i=1:k
return

