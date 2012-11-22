function [x, history] = admm_gl_fs(varargin)
% group_lasso_feat_split  Solve group lasso problem via ADMM feature splitting
%
% [x, history] = group_lasso_feat_split(A, b, p, lambda, rho, alpha);
% 
% solves the following problem via ADMM:
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda sum(norm(x_i))
%
% The input p is a K-element vector giving the block sizes n_i, so that x_i
% is in R^{n_i}.
% 
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and 
% dual residual norms, and the tolerances for the primal and dual residual 
% norms at each iteration.
% 
% rho is the augmented Lagrangian parameter. 
%
% alpha is the over-relaxation parameter (typical values for alpha are 
% between 1.0 and 1.8).
%
% This version is a (serially) distributed, feature splitting example.
%
% Author: Based on example code by Boyd et al [1-2]. 
%         Tim Mullen, 2012, SCCN/INC/UCSD
%
% References:
% [1] http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
% [2] http://www.stanford.edu/~boyd/papers/admm/group_lasso/group_lasso_feat_split.html

g = arg_define(3,varargin, ...
                arg_norep({'A','DesignMatrix'},mandatory,[],'The design matrix. This is the data matrix.'), ...
                arg_norep({'b','TargetVector'},mandatory,[],'The target vector'), ...
                arg_norep({'blks','Blocks'},mandatory,[],'Array containing the lengths of each block'), ...
                arg({'lambda','RegularizationParam'},[],[0 Inf],'Regularization parameter (lambda)'), ...
                arg({'rho','AugmentedLagrangianParam'},1.0,[0 Inf],'Initial value of augmented Lagrangian parameter (rho)'), ...
                arg({'alpha','OverRelaxationParam','Alpha'},1.7,[0 Inf],'Over-relaxation parameter (alpha). Typical values are between 1.5 and 1.8'), ...
                arg({'z_init','InitialState','initAR'},[],[],'Initial state vector','shape','matrix'), ...
                arg({'verb','Verbosity'},false,[],'Verbose output'), ...
                arg({'max_iter','MaxIterations'},1000,[],'Maximum number of iterations'), ...
                arg({'compute_objval','ComputeObjectiveValue'},false,[],'Compute objective value function. Slower processing. Useful mainly for diagnostics'), ...
                arg({'rho_update','RhoUpdate'},false,[],'Update Rho. Whether to update rho dynamically according to 3.4.1 in [1]. Note, this can sometimes cause r_norm, s_norm to "blow up"'), ...
                arg({'rho_cutoff','RhoUpdateThreshold'},10.0,[],'Rho update threshold.'), ...
                arg({'rho_incr','RhoUpdateIncr'},2.0,[],'Rho update increment factor.'), ...
                arg({'rho_decr','RhoUpdateDecr'},2.0,[],'Rho update decrement factor.'), ...
                arg({'lambda_update','LambdaUpdate'},false,[],'Update Lambda. Whether to decrease lambda (regularization) dynamically if convergence is slow.'), ...
                arg({'lambda_update_thresh','LambdaUpdateThreshold'},10^-5,[],'Update lambda if the change in r_norm is less than this','type','denserealdouble'), ...
                arg({'lambda_update_count','LambdaUpdateCount'},10,[],'Number of iterations before updating lambda. Update lambda if r_norm has not significantly changed after this many iterations.'), ...
                arg({'lambda_update_factor','LambdaUpdateFactor'},10,[],'Update factor by which to divide lambda.'), ...
                arg({'abstol','AbsoluteTolerance'},1e-4,[0 Inf],'Primal/Dual absolute tolerance'), ...
                arg({'reltol','RelativeTolerance'},1e-2,[0 Inf],'Primal/Dual relative tolerance'), ...
                arg_sub({'x_update_opts','X_UpdateSolverOptions'},{},...
                {...
                arg({'max_iters','MaxIterations'},100,[],'Number of iterations.'), ...
                arg({'tol','Tolerance'},1e-6,[],'Stopping tolerance. Default value is 1e-6. However, this may be relaxed substantially as described in section 4.3 of [1].'), ...
                },'Options for iterative bisection solver in x-update. Here x is the unknown parameter vector') ...
                );
            
                
arg_toworkspace(g);

if verb
    t_start = tic;
end

%% Data initialization
[m, n] = size(A);

% check that blks divides in to n
if (rem(n,blks) ~= 0)
    error('invalid block size');
end
% number of subsystems
N = n/blks;

x = zeros(blks,N);
u = zeros(m,1);
if isempty(z_init)
    z = zeros(m,1);  
else
    z = z_init;
end
Axbar = zeros(m,1);

zs = zeros(m,N);
Aixi = zeros(m,N);

% extract blocks for fast parallel computation
Ai = cell(1,N);
for i=1:N
    Ai{i} = A(:,(i-1)*blks+1:i*blks);
end

if verb
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

history = struct(...
            'r_norm'  ,   nan(max_iter,1), ...
            's_norm'  ,   nan(max_iter,1), ...
            'eps_pri' ,   nan(max_iter,1), ...
            'eps_dual',   nan(max_iter,1), ...
            'objval'  ,   nan(max_iter,1), ...
            'xup_iter',   nan(max_iter,1)  ...
            );
        
%% ADMM solver
lambda_counter = 0;

% pre-factor (in parallel)
[V D] = deal(cell(1,N));
parfor i = 1:N
    [Vi,Di] = eig(Ai{i}'*Ai{i});
    V{i} = Vi;
    D{i} = diag(Di);
    
    % in Matlab, transposing costs space and flops
    % so we save a transpose operation everytime
    At{i} = Ai{i}';
end   

for k = 1:max_iter
    
    % x-update (to be done in parallel)
    LAMBDA = lambda;
    RHO    = rho;
    NITER  = x_update_opts.niter;
    TOL    = x_update_opts.tol;
    parfor i = 1:N,
%         Ai = A(:,(i-1)*blks + 1:i*blks);
        xx = x_update(Ai{i}, At{i}, Aixi(:,i) + z - Axbar - u, ...
                      LAMBDA/RHO, V{i}, D{i}, NITER, TOL);
        x(:,i) = xx;
        Aixi(:,i) = Ai{i}*x(:,i);     
    end
    
    % z-update
    zold = z;
    Axbar = 1/N*A*vec(x);
    
    Axbar_hat = alpha*Axbar + (1-alpha)*zold;
    z = (b + rho*(Axbar_hat + u))/(N+rho);
    
    % u-update
    u = u + Axbar_hat - z;
    
    % compute the dual residual norm square
    s = 0; q = 0;
    zsold = zs;
    zs = z*ones(1,N) + Aixi - Axbar*ones(1,N);
    for i = 1:N,
        % dual residual norm square
        s = s + norm(-rho*At{i}*(zs(:,i) - zsold(:,i)))^2;
        % dual residual epsilon
        q = q + norm(rho*At{i}*u)^2;
    end

    % diagnostics, reporting, termination checks
    if g.compute_objval
        history.objval(k)  = objective(A, b, lambda, N, x, z);
    end
    history.r_norm(k)  = sqrt(N)*norm(z - Axbar);
    history.s_norm(k)  = sqrt(s);
    
    history.eps_pri(k) = sqrt(n)*abstol + reltol*max(norm(Aixi,'fro'), norm(-zs, 'fro'));
    history.eps_dual(k)= sqrt(n)*abstol + reltol*sqrt(q);

    
    if verb
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end
    
    if history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k);
        break
    end

    % update rho              
    if rho_update
        if history.r_norm(k) > rho_cutoff * history.s_norm(k)
            rho = rho .* rho_incr;
            u   = u ./ rho_incr;
            if verb, fprintf('  incrementing rho to %5.5g.\n',rho); end
        elseif history.s_norm(k) > rho_cutoff * history.r_norm(k)
            rho = rho ./ rho_incr;
            u   = u .* rho_incr;
            if verb, fprintf('  decrementing rho to %5.5g.\n',rho); end
        end
    end
    
    if lambda_update && k > 1
        % check if convergence has not improved in a while
        % and, if so, relax shrinkage a bit
        if  abs(history.r_norm(k)-history.r_norm(k-1)) < lambda_update_thresh ...
            && abs(history.s_norm(k)-history.s_norm(k-1)) < lambda_update_thresh
            lambda_counter = lambda_counter + 1;
            
            if lambda_counter > lambda_update_count
                lambda = lambda/lambda_update_factor;
                lambda_counter = 0;
                if verb, fprintf('  relaxing lambda to %0.7g.\n',lambda); end
            end
        end
    end
    
    if verb && mod(k,15)
    	fprintf('\n%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n\n',  ...
                'iter', 'r norm', 'eps pri', 's norm',      ...
                'eps dual', 'objective');
    end
    
end % loop over k=1:max_iter

if k>=max_iter && verb
    fprintf('glADMM:Maximum number of iterations (%d) reached\n',max_iter);
end

if verb
    toc(t_start);
end

end


%% objective function
function p = objective(A, b, lambda, N, x, z)
    p = ( 1/2*sum_square(N*z - b) + lambda*sum(norms(x)) );
end

%% x-update function
function [x i] = x_update(A, At, b, kappa, V, D, niter, tol)
    [m,n] = size(A);

    q = At*b;

    if (norm(q) <= kappa)
       x = zeros(n,1);
    else
        % bisection on t
        lower = 0; upper = 1e10;
        for i = 1:niter,
            t = (upper + lower)/2;

            x = V*((V'*q)./(D + t));
            if t > kappa/norm(x),
                upper = t;
            else
                lower = t;
            end
            if (upper - lower <= tol)
                break;
            end
        end
    end
end

