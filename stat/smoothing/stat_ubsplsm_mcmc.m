function [FIT MCMC_LastState]=stat_ubsplsm_mcmc(varargin)
% Y is a cell array where Y{i} is the T x 1 vector of connectivities for
% the ith channel pair of the jth subject
% Y = {s1(1,1) s1(1,2) s1(1,3) ... s1(2,1) s1(2,2) s1(2,3) ... s2(1,1) s2(1,2) ...}
% If there are multiple subjects, then each subject's channel pairs are 
% simply appended as additional cells of Y{:}. Thus for NC channels and NS
% subjects, Y is at most of dimension [1 x NS*NC^2] (if diagonals are included)

% FIT contains the MCMC estimates
% MCMC_LastState contains the final state of the Gibbs sampler

if ~exist('verb','var')
    verb = 2;
end

g = arg_define([0 1],varargin, ...
    arg_norep({'Y','TSData'},mandatory,[],sprintf(['Cell array of data to smooth.\n' ...
              'Generally, Y is a cell array where Y{i,s} is the T x 1 vector of time-varying (or freq-varying) connectivity for the ith channel pair of the sth subject.\n' ...
              'e.g. Y = {s1(1,1) s1(1,2) s1(1,3) ... s1(N,1) s1(N,2) s1(N,3) ... \n' ...
              '          s2(1,1) s2(1,2) s2(1,3) ... } \n'])), ...
    arg({'K','Knots'},5,[],'Positions of spline knots along frequency dimension (Hz). If K is a scalar, then K knots are evenly spaced from first to last frequency. A good heuristic is one knot every 5%','shape','row'), ...
    arg({'Q','FPCABasisDim','fpcaBasisDim'},4,[0 Inf],'Number of FPCA basis functions.'), ...
    arg({'smoothingLayout','MatrixElementsToSmooth'},{'diagonals','off-diagonals'},{'off-diagonals'},'Which parts of the matrix to smooth. Diagonals (e.g. auto-connectivity) and off-diagonals (e.g. cross-connectivity) will be smoothed separately','type','logical'), ...
    arg({'niters','nMCMCiters','NumMcmcIters','niter'},1000, [1 Inf], 'Number of MCMC iterations for spline fitting'), ...
    arg({'basisCoeffVarPrior'},1000,[eps Inf],'Variance of basis coefficient gaussian prior. Larger --> more wiggling allowed'), ...
    arg({'noiseVarPriorShape'},0.01,[eps Inf],'Shape (D.O.F) of noise variance prior. This is the "alpha" parameter of the inverse gamma prior distribution. Increasing noiseVarPriorShape --> decreased variance of noise variance distribution.'), ...
    arg({'noiseVarPriorScale'},0.01,[eps Inf],'Scale parameter of noise variance prior. This is the "theta" (1/beta) parameter of inverse gamma prior distribution. Increasing noiseVarPriorScale --> right-shift of distribution --> (increase in expected noise variance). In general MEAN(noiseVariance) = noiseVarPriorScale/noiseVarPriorShape and MODE(noiseVariance) = noiseVarPriorScale/(noiseVarPriorShape-1) for noiseVarPriorShape>=1.'), ...
    arg({'initNoiseVariance'},0.1,[eps Inf],'Initial noise variance'), ...
    arg_norep({'MCMC_InitState'},struct([]),[],'Object containing initial state of Gibbs sampler'));
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...

arg_toworkspace(g);

% make Y a column vector
Y = Y(:);

R=size(Y,1);        % number of time-series to smooth across all subjects
T=size(Y{1},1);     % number of time points

%% Determine knot locations
% -------------------------------------------------------------------------
if isscalar(K)
    idx   = 1:T;
    order = 4;
    K = K+order-1;
    knots = quantile(idx,(0:1:(K-order+1)) / (K-order+1));
    knots = knots(:)';
else
    % knots are already provided
    knots = K(:);
    K = length(knots)+2;
end
      

%% Initialize Parameters
% -------------------------------------------------------------------------
ALPHA        = zeros(Q,R,niters+1);      % For each channel, the Q-dimensional vector of spline regression weights
ALPHA_BAR    = zeros(Q,niters+1);        % Multivariate gaussian prior distribution for ALPHA
SIGMA_EPS    = ones(1,niters+1);         % noise variance
SIGMA_EPS(1) = initNoiseVariance;        % initial noise variance
FIT          = zeros(T,R,niters+1);      % Smoothed measure estimates
THETA        = zeros(K,Q,niters+1);      % Q FPCA component vectors
        

if ~isempty(MCMC_InitState)
    
    % set initial state to user-provided
    ALPHA(:,:,1)        = MCMC_InitState.ALPHA(:,:,end);
    ALPHA_BAR(:,:,1)    = MCMC_InitState.ALPHA_BAR(:,:,end);
    SIGMA_EPS(1)        = MCMC_InitState.SIGMA_EPS(end);
    THETA(:,:,1)        = MCMC_InitState.THETA(:,:,end);
    phi_t               = MCMC_InitState.phi_t;
    
else
    
    % Construct orthonormal basis functions
    % -------------------------------------------------------------------------
    phi_t = stat_ubsplsm_mkspl(knots,T,4,verb);

    % Initialize FPCA
    % -------------------------------------------------------------------------
    if verb==2
        multiWaitbar('Initializing FPCA','Reset','Color',hlp_getNextUniqueColor);
    end

    % Initialize the FPCA components to random, orthonormal vectors
    THETA(:,1,1) = mvnrnd(zeros(K,1),eye(K));
    THETA(:,1,1) = THETA(:,1,1)/norm(THETA(:,1,1));

    for q = 2:Q                                                                     % [!] should be able to replace this whole section with a single line of orth(randn(...)). We only need orthonormal random vectors
        if verb==2
            multiWaitbar('Initializing FPCA',q/Q);
        end
        THETA(:,q,1) = mvnrnd(zeros(K,1),eye(K));

        % orthonormalization of random vector
        for r = 1:(q-1)

            THETA(:,q,1) = THETA(:,q,1) ...
                           -(THETA(:,q,1)'*THETA(:,q-r,1)/norm(THETA(:,q-r,1).^2)) ...
                           * THETA(:,q-r,1);
        end
        THETA(:,q,1) = THETA(:,q,1)/norm(THETA(:,q,1));
    end
    
end

% compute temporal sum of splice basis function covariance matrices:
% e.g. Sigma_phi_sum := sum_t { phi_t(:,t)*phi_t(:,t)' } = phi_t*phi_t';
Sigma_phi_sum = phi_t*phi_t';

%% run MCMC to smooth data
% -------------------------------------------------------------------------
if verb==2
    multiWaitbar('Gibbs sampling','Reset','Color',hlp_getNextUniqueColor);
end

% constants
% -------------------------------------------------------------------------
Iq              = eye(Q);
Iq_sp           = speye(Q);
Iqk             = eye(Q*K);
Sigma_alpha_bar = Iq/(R + 1/basisCoeffVarPrior);

for iter=1:niters
   

    % Draw ALPHA (spline regression coefficients (weights))
    % ---------------------------------------------------------------------
    
    % compute spline regression coeffs inverse covariance matrix
    Sigma_alpha_i     = Iq + (THETA(:,:,iter)'*Sigma_phi_sum*THETA(:,:,iter))/SIGMA_EPS(iter);
    Sigma_alpha_i     = covfixer(Sigma_alpha_i);
    
    % invert inverse covariance matrix to produce cov mat
    Sigma_alpha_i     = inverse(Sigma_alpha_i);
    Sigma_alpha_i_dbl = double(Sigma_alpha_i);
    
    % pre-compute product
    tmpprod = THETA(:,:,iter)'*phi_t/SIGMA_EPS(iter);                      
    for pair=1:R  % for each channel pair
        
        % compute regression coeffs mean
        mu_alpha_i = Sigma_alpha_i*(ALPHA_BAR(:,iter) + tmpprod*Y{pair});
        
        % draw spline regression coefficients for this channel pair
        ALPHA(:,pair,iter+1) = mvnrnd(mu_alpha_i,Sigma_alpha_i_dbl)';
    end
    
    
    % Draw ALPHA_BAR (hyperparameter for ALPHA gaussian prior mean)
    % ---------------------------------------------------------------------
    mu_alpha_bar        = sum(ALPHA(:,:,iter+1),2);
    mu_alpha_bar        = Sigma_alpha_bar*mu_alpha_bar;
    ALPHA_BAR(:,iter+1) = mvnrnd(mu_alpha_bar,Sigma_alpha_bar)';
    
    % Draw SIGMA_EPS (noise/residual variance)
    % ---------------------------------------------------------------------
    pi_eps_shape = noiseVarPriorShape + R*T/2;       % prior distribution shape 
    
    % compute the sum-squared error of residuals (data-fit)^2
    pi_eps_scale = noiseVarPriorScale;               % prior distr. scale param
    tmpprod      = phi_t'*THETA(:,:,iter);      % precompute product        
    for pair=1:R
        pi_eps_scale = pi_eps_scale ...
               + sum((Y{pair}-tmpprod*ALPHA(:,pair,iter+1)).^2 / 2);            
    end
    % draw noise variance estimates from inverse gamma
    % note that we assume the noise variance is identical for all pairs         
    SIGMA_EPS(:,iter+1)=1/gamrnd(pi_eps_shape,1/pi_eps_scale);                  
    
    % Draw THETA (FPCA component vectors)
    % ---------------------------------------------------------------------
    
    % NOTE: the Sigma_theta on the inner loop is updated for continuing 
    % iterations of the outer loop...in other words, if there are multiple
    % subjects/channel pairs, the final Sigma_theta value depends on all
    % subjects/channel pairs
    
    % initialize mean and covariance parameters
    Sigma_theta = diag(1/basisCoeffVarPrior); % Iqk/basisCoeffVarPrior;
    mu_theta    = zeros(Q*K,1);
    
    for i = 1:R  % for each channel pair
        Alpha_i    = ALPHA(:,i,iter+1)';
        Alpha_i_sq = Alpha_i'*Alpha_i;
        for t = 1:T  % for each time window                                     % [!] we should be able to get rid of this loop
            Phi_D_t     = kron(Iq_sp,phi_t(:,t)');                              % [!] this can be pre-computed and moved to outside the loop
            Sigma_theta = Sigma_theta + Phi_D_t'*Alpha_i_sq*Phi_D_t;
            mu_theta    = mu_theta + Phi_D_t'*Alpha_i'*Y{i}(t);
        end
    end
    
    Sigma_theta = Sigma_theta/SIGMA_EPS(iter+1);

    % DEBUG: PROFILING
    
    % enforce valid covariance matrix before and after inversion
    Sigma_theta = covfixer(Sigma_theta);
    Sigma_theta = Iqk/Sigma_theta;
    Sigma_theta = covfixer(Sigma_theta);
    
    % compute mu (mean)
    mu_theta    = mu_theta/SIGMA_EPS(iter+1);
    mu_theta    = Sigma_theta*mu_theta;
    
    % draw theta
    THETA(:,:,iter+1) = mvnrnd(mu_theta,Sigma_theta)';
    THETA(:,:,iter+1) = reshape(THETA(:,:,iter+1),K,Q);
    
    % Calculate smoothed fits
    % ---------------------------------------------------------------------
    FIT(:,:,iter+1) = phi_t'*THETA(:,:,iter+1)*ALPHA(:,:,iter+1);
    
    if verb==2
        multiWaitbar('Gibbs sampling',iter/niters);
    end
    
end


if nargout>2
    % store final state of Gibbs sampler
    MCMC_LastState = struct( ...
     'THETA'            , THETA(end),       ...   % FPCA component vectors
     'ALPHA'            , ALPHA(end),       ...   % Spline regression weights
     'ALPHA_BAR'        , ALPHA_BAR(end),   ...   % Gaussian prior for ALPHA
     'SIGMA_EPS'        , SIGMA_EPS(end),   ...   % Noise variance
     'phi_t'            , phi_t,            ...   % Spline basis functions
     'Sigma_theta'      , Sigma_theta,      ...   % THETA cov mat (mvnrnd)
     'mu_theta'         , mu_theta,         ...   % THETA mean (mvnrnd)
     'Sigma_alpha_i'    , Sigma_alpha_i_dbl,...   % ALPHA cov mat (mnvnrnd)
     'mu_alpha_i'       , mu_alpha_i,       ...   % ALPHA mean (mnvnrnd)
     'Sigma_alpha_bar'  , Sigma_alpha_bar,  ...   % ALPHA_BAR cov mat (mnvnrnd)
     'mu_alpha_bar'     , mu_alpha_bar      ...   % ALPHA_BAR mean (mnvnrnd)
     );
end
                     

if verb==2
    % cleanup waitbars
    multiWaitbar('Gibbs sampling'   , 'Close');
    multiWaitbar('Initializing FPCA', 'Close');
    multiWaitbar('Building Splines' , 'Close');
end





% DEBUG SECTION

% -------------------------------------------------------------------------
% PROFILING
% -------------------------------------------------------------------------
%     fprintf('\n')
%
%
%
%     Alpha_i_sq_sp = sparse(Alpha_i_sq);
%
%
%
%     T = 500;
%
%     tic
%     for k=1:T
%         Phi_D_t  =kron(Iq,phi_t(:,1)');
%         Phi_D_t'*Alpha_i_sq*Phi_D_t;
%     end
%     fprintf('full:\t\t%0.5g\n', toc/T)
%
%     tic
%     for k=1:T
%         Phi_D_t  =kron(Iq,phi_t(:,1)');
%         Phi_D_t_tr  =kron(phi_t(:,1),Iq);
%         Phi_D_t_tr*Alpha_i_sq*Phi_D_t;
%     end
%     fprintf('full (no T):\t\t%0.5g\n', toc/T)
%
%     tic
%     for k=1:T
%         Phi_D_sp =kron(Iq_sp,phi_t(:,1)');
%         Phi_D_sp'*Alpha_i_sq*Phi_D_sp;
%     end
%     fprintf('sparse:\t\t%0.5g\n', toc/T)
%
%     tic
%     for k=1:T
%         Phi_D_sp =kron(Iq_sp,phi_t(:,1)');
%         Phi_D_sp_tr =kron(phi_t(:,1),Iq_sp);
%         Phi_D_sp_tr*Alpha_i_sq*Phi_D_sp;
%     end
%     fprintf('sparse (no T):\t\t%0.5g\n', toc/T)
%
%
%     % fast version
% -------------------------------------------------------------------------
%     Sigma_theta2=eye(Q*K)/basisCoeffVarPrior;
%     mu_theta2=zeros(Q*K,1)
%
%     QK = Q*K;
%     Phi_D_t = zeros(QK*T,Q);
%     for t=1:T
%         % form kroneker matrix (can perhaps be optimized)
%         Phi_D_t((t-1)*QK+1:QK*((t-1)+1),:) = kron(Iq,phi_t(:,t)');
%     end
%
%     for(i=1:R)
%         Alpha_i2=ALPHA(:,i,iter+1)';
%         Alpha_i_sq = Alpha_i2'*Alpha_i2;
%
%         tmp = kron(Phi_D_t'*Alpha_i_sq,Phi_D_t);
%
%         Sigma_theta2=Sigma_theta2+tmp/SIGMA_EPS(iter+1);
%
% %         Sigma_theta = sum(reshape(Sigma_theta,[Q,
% %         mu_theta=mu_theta+Phi_D_t'*Alpha_i'*Y{i}(t)/SIGMA_EPS(iter+1);
%     end
