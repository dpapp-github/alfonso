function [delta, probData] = linsolveB(soln, probData, RHS)
% This method implements a block linear solver for the Newton systems.
% --------------------------------------------------------------------------
% USAGE of "linsolveB"
% delta = linsolveB(soln, probData, RHS)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% RHS:          right-hand side of the Newton system
%
% OUTPUT
% delta:	    computed Newton direction
% --------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Version: 2024/07/15
% 
% This code has been developed and tested with Matlab R2023b.
% -------------------------------------------------------------------------
       
    A  = probData.A;
    At = probData.At;
    b  = probData.b;
    c  = probData.c;
    [m, n] = size(A);

    mu     = soln.mu;
    tau    = soln.tau;
    Hi     = soln.Hi;
    AHiAti = computeAHiAti(soln.Li, probData);
    %?%AHiAti = soln.AHiAti;
    
    ry     = RHS(1:m);
    rx     = RHS(m+(1:n));
    rtau   = RHS(m+n+1);
    rs     = RHS(m+n+1+(1:n));
    rkappa = RHS(end);

    % Schur complement solver for the system
    % [0   A  -b]  [dy]       [ry]
    % [-A' H   c]  [dx]    =  [rx+rs]
    % [b' -c'  h]  [dtau]     [rtau+rkappa]

    % equivalently
    % [ H  -A'  c]  [dx]       [rx+rs]
    % [-A   0   b]  [dy]    =  [-ry]
    % [ c' -b' -h]  [dtau]     [-rtau-rkappa]

    % The next few steps: dx, dy, dtau are modified step-by-step.
    % Notation for the comments below: inv(LHS) = UDL block factorization
    dx   = rx+rs;
    dy   = -ry;
    dtau = -rtau-rkappa;
    
    % Precompute dense rows/columns of L and U.
    % Slightly wasteful for better modularity. (Hi(c) is computed twice.)
    [u13, u23] = HAA0sol(c, b, probData, mu, Hi, AHiAti);
    [l31, l32] = HAA0sol(c, -b, probData, mu, Hi, AHiAti);
    %l31 = u13; not true, although l32 = -u23; ?

    alpha = -mu/tau^2 - c'*u13 + b'*u23;

    % Multiply by L.
    dtau = dtau - l31'*dx - l32'*dy;

    % Multiply by D.
    [dx, dy] = HAA0sol(dx, dy, probData, mu, Hi, AHiAti);
    dtau = dtau/alpha;

    % Multiply by U.
    dx = dx - u13*dtau;
    dy = dy - u23*dtau;

    delta               = zeros(m+2*n+2, 1);
    delta(1:m)          = dy;
    delta(m+n+1)        = dtau;
    delta(m+(1:n))      = dx;
    delta(m+n+1+(1:n))  = -rx - At*dy + dtau*c;   % ds
    delta(end)          = -rtau + b'*dy - c'*dx;  % dkappa
    
return

% Computes [H -A'; -A 0] \ [u;v]
function [u,v] = HAA0sol(u, v, probData, mu, Hi, AHiAti)

    % Notation for the comments below: inv([H -A'; -A 0]) = UDL block factorization

    % Multiply by L.
    if ~iscell(Hi)
        Hiu = Hi(u)/mu;
    else
        Kdims = probData.Kdims;

        idx = 0;
        Hiu = zeros(length(u),1);
        for k=1:length(Kdims)
            Hiu(idx+1:idx+Kdims(k)) = Hi{k}(u(idx+1:idx+Kdims(k)))/mu;
            idx = idx+Kdims(k);
        end
    end
    v = v + probData.A*Hiu;

    % Multiply by D.
    u = Hiu;
    v = -AHiAti(v)*mu;

    % Multiply by U.
    if ~iscell(Hi)
        u = u + Hi(probData.At*v)/mu;
    else
        idx = 0;
        for k=1:length(Kdims)
            u(idx+1:idx+Kdims(k)) = u(idx+1:idx+Kdims(k)) + Hi{k}(probData.As{k}*v(probData.rB(:,k)))/mu;
            idx = idx+Kdims(k);
        end
    end
    %v = v*mu;

return
