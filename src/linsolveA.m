function [delta, probData] = linsolveA(soln, probData, RHS)
% This method implements a block linear solver for the Newton systems.
% --------------------------------------------------------------------------
% USAGE of "linSolveA"
% delta = linSolveA(soln, probData, RHS)
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
    
    if issparse(A)
        b = sparse(b);
        c = sparse(c);
        RHS = sparse(RHS);
    end
    
    mu     = soln.mu;
    tau    = soln.tau;
    Hi     = soln.Hi;
    Li     = soln.Li;
    %AHiAti = soln.AHiAti;
    AHiAt  = computeAHiAt(Li, probData);
    
    ry     = RHS(1:m);
    rx     = RHS(m+(1:n));
    rtau   = RHS(m+n+1);
    rs     = RHS(m+n+1+(1:n));
    rkappa = RHS(end);

    if ~iscell(Hi)
        %LiAt    = Li(At);
        Lic     = Li(c);
        Hic     = Hi(c);
        Hirxrs  = Hi(rx+rs);
    else
        Kdims = probData.Kdims;

        idx = 0;
        Lic    = zeros(n,1);
        Hic    = zeros(n,1);
        Hirxrs = zeros(n,1);

        z = rx+rs;
        for k=1:length(Kdims)
            ck = c(idx+1:idx+Kdims(k));
            zk = z(idx+1:idx+Kdims(k));
            Lic(idx+1:idx+Kdims(k)) = Li{k}(ck);
            Hic(idx+1:idx+Kdims(k)) = Hi{k}(ck);
            Hirxrs(idx+1:idx+Kdims(k)) = Hi{k}(zk);
            idx = idx+Kdims(k);
        end
    end
    AHic    = A*Hic;
    
    RHSdydtau   = [ry; rtau+rkappa] - [A; -c']*Hirxrs/soln.mu;
    
    x = RHSdydtau(1:end-1);
    y = RHSdydtau(end);

    dinv = (mu/tau^2 + Lic.'*Lic/mu)^(-1);
    %dy   = (LiAt'*LiAt/mu + (b+AHic/mu)*dinv*(b-AHic/mu)') \ (x+(b+AHic/mu)*dinv*y);
    dy   = (AHiAt/mu + (b+AHic/mu)*dinv*(b-AHic/mu)') \ (x+(b+AHic/mu)*dinv*y);
    dtau = dinv*((-b + AHic/mu)'*dy + y);
    
    %dx = (Hirxrs + L'\([LiAt, -Lic]*[dy; dtau]))/soln.mu;
    if ~iscell(Hi)
        HiAtdy = Hi(At*dy);
    else
        idx = 0;
        HiAtdy = zeros(n,1);

        z = At*dy;
        for k=1:length(Kdims)
            zk = z(idx+1:idx+Kdims(k));
            HiAtdy(idx+1:idx+Kdims(k)) = Hi{k}(zk);
            idx = idx+Kdims(k);
        end
    end
    dx = (Hirxrs + HiAtdy - Hic*dtau)/soln.mu;
    
    delta               = zeros(m+2*n+2, 1);
    delta(1:m)          = dy;
    delta(m+n+1)        = dtau;
    delta(m+(1:n))      = dx;
    delta(m+n+1+(1:n))  = -rx - [At, -c]*[dy; dtau];
    delta(end)          = -rtau + b'*dy - c'*dx;
    
return
