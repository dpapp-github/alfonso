function e_design()
% Numerical experiments with E-optimal design problems comparing alfonso,
% Mosek 9.2.6 and SCS 2.1.7
%
% -------------------------------------------------------------------------
% EXTERNAL PACKAGES CALLED IN THIS FILE: YALMIP, Mosek, SCS
% Because of these dependencies, this example will NOT run in Octave.

diary('e_design_results.txt');
datetime()
rng(2020,'twister');

for i=1:5 % warning: i=6:10 is too slow for SCS
    
    n = 50*i;
    p = 2*n;

    pars = struct('v', randn(n,p));

    probData = struct('c',[-1;zeros(p,1)],'A',[0,ones(1,p)],'b',1);
    x0 = [0; ones(p,1)/p];
    opts = struct('optimTol',1e-8);
    
    % ALFONSO
    res=alfonso(probData,x0,@ex_oracle_design,pars,opts);
    % print # of iterations and solver time
    fprintf('alf: %d %d %.2f\n\n', n, res.nIterations, res.time);

    % YALMIP PREP FOR MOSEK AND SCS
    x = sdpvar(p,1);
    t = sdpvar(1);
    cons = [sum(x)==1, x>=0, -t*eye(n)+pars.v*diag(x)*pars.v' >= 0];

    % MOSEK
    opts = sdpsettings('solver','mosek','savesolveroutput',1);
    res =  optimize(cons,-t,opts); % maximize t
    % print only solver (not yalmip) time
    fprintf('mos: %d %d %.2f\n\n', n, res.solveroutput.res.info.MSK_IINF_INTPNT_ITER, res.solvertime);

    % SCS
    opts = sdpsettings('solver','scs','savesolveroutput',1,'scs.max_iters',100000,'scs.eps',1e-3);
    res =  optimize(cons,-t,opts); % maximize t
    % print only solver (not yalmip) time
    fprintf('scs: %d %d %.2f\n\n', n, res.solveroutput.info.iter, res.solvertime); 
end

disp('Done!');
diary('off');

return


function [in, g, H, L] = ex_oracle_design(tx, pars)
% This function implements a membership and barrier function oracle for
% the E-optimal design example e_design.m
%
% INPUT
% tx:                   column vector representing [t; x(1); ...; x(n)]
% pars:                 structure with a single field pars.v
%                       pars.v is a two-dimensional array whose ith column
%                       v(:,i) represents the ith design vector (i=1,...,p)

t = tx(1);
x = tx(2:end);
[n,p] = size(pars.v);

% in the cone?
if any(x <= 0)
    in = false;
    g  = NaN;
    H  = NaN;
    L  = NaN;
    return
end

Ax = -t*eye(n) + pars.v*diag(x)*pars.v';

[L,err] = chol(Ax,'lower');
if err > 0
    in = false;
    g  = NaN;
    H  = NaN;
    L  = NaN;
    return
else
    in = true;
end

% tx is in the cone

% compute g and H if required
if nargout > 1
    
    g = [0; -1./x];

    Li = inv(L);
    g(1) = Li(:)'*Li(:);
    w = L\pars.v;
    for i=1:p
       g(i+1) = -w(:,i)'*w(:,i);
    end
    
    % compute H and L if required
    if nargout > 2
        H = diag([0; x.^(-2)]);
        
        invAx = Li'*Li;
        H(1,1) = H(1,1) + invAx(:)'*invAx(:);
        Lws = L' \ w;
        for i=2:p+1
            H(i,1) = H(i,1) - Lws(:,i-1)'*Lws(:,i-1);
        end
        H(1,2:p+1) = H(2:p+1,1)';
        H(2:end,2:end) = H(2:end,2:end) + (w'*w).^2;
        
        if nargout > 3
            [L,err] = chol(H,'lower');
            if err > 0
                in = false;
                g  = NaN;
                H  = NaN;
                L  = NaN;
                return
            end
        end
    end
end
    
return
