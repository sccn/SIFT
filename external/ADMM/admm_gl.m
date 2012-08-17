function [z, history] = admm_gl(varargin)  %A, b, lambda, p, rho, alpha, z_init, verb, max_iter

% group_lasso  Solve group lasso problem via ADMM
%
% [x, history] = group_lasso(A, b, blks, ...);
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
%
% More information can be found in the papers found in [1-2]
%
% Based on example code by Boyd et al [1-2].
% Modified by: Tim Mullen, 2012, SCCN/INC/UCSD
%
% References:
% [1] http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
% [2] http://www.stanford.edu/~boyd/papers/admm/group_lasso/group_lasso_example.html

g = arg_define(3,varargin, ...
                arg_norep({'A','DesignMatrix'},mandatory,[],'The design matrix'), ...
                arg_norep({'b','TargetVector'},mandatory,[],'The target vector'), ...
                arg_norep({'blks','Blocks'},mandatory,[],'Array containing the lengths of each block'), ...
                arg({'lambda','RegularizationParam'},[],[0 Inf],'Regularization parameter'), ...
                arg({'rho','AugmentedLagrangianParam'},1.0,[0 Inf],'Initial value of augmented Lagrangian parameter'), ...
                arg({'Alpha','OverRelaxationParam','alpha'},1.7,[0 Inf],'Over-relaxation parameter. Typical values are between 1.5 and 1.8'), ...
                arg({'z_init','InitialState','initAR'},[],[],'Initial state vector','shape','matrix'), ...
                arg({'verb','Verbosity'},false,[],'Verbose output'), ...
                arg({'max_iter','MaxIterations'},1000,[],'Maximum number of iterations'), ...
                arg({'compute_objval','ComputeObjectiveValue'},false,[],'Compute objective value function. Slower processing. Useful mainly for diagnostics'), ...
                arg({'rho_update','RhoUpdate'},false,[],'Update Rho. Whether to update rho dynamically according to 3.4.1 in [1]. Note, this can sometimes cause r_norm, s_norm to "blow up"'), ...
                arg({'rho_cutoff','RhoUpdateThreshold'},10.0,[],'Rho update threshold.'), ...
                arg({'rho_incr','RhoUpdateIncr'},2.0,[],'Rho update increment factor.'), ...
                arg({'rho_decr','RhoUpdateDecr'},2.0,[],'Rho update decrement factor.'), ...
                arg({'lambda_update','LambdaUpdate'},false,[],'Update Lambda. Whether to decrease lambda (regularization) dynamically if convergence is slow.'), ...
                arg({'lambda_update_thresh','LambdaUpdateThreshold'},10^-5,[],'Update lambda if the change in r_norm is less than this'), ...
                arg({'lambda_update_count','LambdaUpdateCount'},10,[],'Number of iterations before updating lambda. Update lambda if r_norm has not significantly changed after this many iterations.'), ...
                arg({'lambda_update_factor','LambdaUpdateFactor'},10,[],'Update factor by which to divide lambda.') ...
                );
                
arg_toworkspace(g);


if verb
    t_start = tic;
end

%% Global constants and defaults
compute_objval = nargout > 1 || compute_objval;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;

%% select lambda using heuristic
% example based on [2]
if isempty(lambda)

    nblks = length(blks);
    
    if g.verb, fprintf('Using heuristic lambda selection...'); end
        
    cum_part = cumsum(blks(1:nblks-1));
    
    % guess regularization param for 
    % each group of coefficients
    K = nblks-1;
    start_ind = 1;
    lambdas = zeros(1,K);
    
    for i = 1:K
        sel = start_ind:cum_part(i);
        lambdas(i) = norm(A(:,sel)'*b);
        start_ind = cum_part(i) + 1;
    end
    lambda_max = max(lambdas);

    % regularization parameter as fraction of 
    % maximum group regularization parameter
    lambda = 0.1*lambda_max;
    
    if g.verb, fprintf('...lambda set to %0.10g\n',lambda); end
end

%% Data preprocessing

[m, n] = size(A);

% save a matrix-vector multiply
Atb = A'*b;
% check that sum(blks) = total number of elements in x
if (sum(blks) ~= n)
    error('invalid partition');
end

% cumulative partition
cum_part = cumsum(blks);

%% ADMM solver

x = zeros(n,1);
u = zeros(n,1);
if isempty(z_init)
    z = zeros(n,1);  
else
    z = z_init;
end

history = struct(...
            'r_norm'  ,   nan(max_iter,1), ...
            's_norm'  ,   nan(max_iter,1), ...
            'eps_pri' ,   nan(max_iter,1), ...
            'eps_dual',   nan(max_iter,1), ...
            'objval'  ,   nan(max_iter,1)  ...
            );
            
% pre-factor
[L U] = factor(A, rho);

if verb
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

lambda_counter = 0;

for k = 1:max_iter

    % x-update
    q = Atb + rho*(z - u);    % temporary value
    if( m >= n )    % if skinny
       x = U \ (L \ q);
    else            % if fat
       x = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
    end

    % z-update
    zold = z;
    start_ind = 1;
    x_hat = Alpha*x + (1-Alpha)*zold;
    for i = 1:length(blks),
        sel = start_ind:cum_part(i);
        z(sel) = shrinkage(x_hat(sel) + u(sel), lambda/rho);
        start_ind = cum_part(i) + 1;
    end
          
    % diagnostics, reporting, termination checks
    history.r_norm(k)  = norm(x - z);
    history.s_norm(k)  = norm(-rho*(z - zold));
    
    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);

    if compute_objval
        history.objval(k)  = objective(A, b, lambda, cum_part, x, z);
    end
    
    if verb
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end

    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
         break;
    end

    % update rho              
    if rho_update
        if history.r_norm(k) > rho_cutoff * history.s_norm(k)
            if verb, disp('  incrementing rho.'); end
            rho = rho .* rho_incr;
%             u = u ./ rho_incr;
        elseif history.s_norm(k) > rho_cutoff * history.r_norm(k)
            if verb, disp('  decrementing rho.'); end
            rho = rho ./ rho_incr;
%             u = u .* rho_incr;
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
                if verb, fprintf('   relaxing lambda to %0.7g.\n',lambda); end
            end
            
        end
    end
    
end

if k>=max_iter && verb
    fprintf('glADMM:Maximum number of iterations (%d) reached\n',max_iter);
end

if verb
    toc(t_start);
end

end

function p = objective(A, b, lambda, cum_part, x, z)
    obj = 0;
    start_ind = 1;
    for i = 1:length(cum_part),
        sel = start_ind:cum_part(i);
        obj = obj + norm(z(sel));
        start_ind = cum_part(i) + 1;
    end
    p = ( 1/2*sum((A*x - b).^2) + lambda*obj );
end

function z = shrinkage(x, kappa)
    z = pos(1 - kappa/norm(x))*x;  % note: need function pos()
end

function [L U] = factor(A, rho)
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A + rho*speye(n), 'lower' );
    else            % if fat
       L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
    end
    
    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end

function q = pos(x)
    q = max(x,0);
end

