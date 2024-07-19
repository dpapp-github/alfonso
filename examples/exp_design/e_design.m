function e_design()
% Numerical experiments with E-optimal design problems comparing alfonso,
% Mosek 9.2.6 and SCS 2.1.7
%
% -------------------------------------------------------------------------
% EXTERNAL PACKAGES CALLED IN THIS FILE: YALMIP, Mosek, SCS
% -------------------------------------------------------------------------

datetime()
rng(2020,'twister');

for i=1:5 % warning: i=6:10 is too slow for SCS
    
    n = 50*i;
    p = 2*n;

    V = randn(n,p);
    
    
    % ALFONSO
    probData = struct('c',[-1;zeros(p,1)],'A',[0,ones(1,p)],'b',1);
    x0 = [0; ones(p,1)/p];
    pars = struct('ext', true, 'nonneg', true, 'Vs', {{V}}, 'ws', ones(p,1));  % cone parameters
    opts = struct('optimTol',1e-8,'linSolveFun',@linsolveB);                   % solver options

    res = alfonso(probData,x0,@gH_rk1LMI,pars,opts);
    % print # of iterations and solver time
    fprintf('alf: %d %d %.2f\n\n', n, res.nIterations, res.time);


    % YALMIP PREP FOR MOSEK AND SCS
    x = sdpvar(p,1);
    t = sdpvar(1);
    cons = [sum(x)==1, x>=0, -t*eye(n)+V*diag(x)*V' >= 0];

    % MOSEK
    opts = sdpsettings('solver','mosek','savesolveroutput',1);
    res =  optimize(cons,-t,opts); % maximize t
    % print only solver (not yalmip) time
    fprintf('mos: %d %d %.2f\n\n', n, res.solveroutput.res.info.MSK_IINF_INTPNT_ITER, res.solvertime);

    %%% SCS
    %%opts = sdpsettings('solver','scs','savesolveroutput',1,'scs.max_iters',100000,'scs.eps',1e-3);
    %%res =  optimize(cons,-t,opts); % maximize t
    %%% print only solver (not yalmip) time
    %%fprintf('scs: %d %d %.2f\n\n', n, res.solveroutput.info.iter, res.solvertime); 
end

disp('Done!');

return

