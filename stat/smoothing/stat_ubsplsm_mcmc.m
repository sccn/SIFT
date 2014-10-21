function [FIT MCMC_LastState] = stat_ubsplsm_mcmc(varargin)
% Y is a cell array where Y{i} is the T x 1 vector of time-varying (or freq-varying) 
% connectivity for the kth channel pair of the sth subject.
% Y = {s1(1,1) s1(1,2) s1(1,3) ... s1(2,1) s1(2,2) s1(2,3) ... s2(1,1) s2(1,2) ...}
% If there are multiple subjects, then each subject's channel pairs are 
% simply appended as additional cells of Y{:}. Thus for NC channels and NS
% subjects, Y is at most of dimension [1 x NS*NC^2] (if diagonals are included)
%
%
% --------------------------------------------------------------------------------------------------------------------------
% Input                  Information                                                                                        
% --------------------------------------------------------------------------------------------------------------------------
% TSData:                Cell array of data to smooth                                                                       
%                        Generally, Y is a cell array where Y{i} is the T x 1 vector of time-varying (or freq-varying)      
%                        connectivity for the kth channel pair of the sk_th subject.                                          
%                        e.g. Y = {s1(1,1) s1(1,2) s1(1,3) ... s1(N,1) s1(N,2) s1(N,3) ...                                  
%                                  s2(1,1) s2(1,2) s2(1,3) ... }                                                            
%                        Y can also be a [1 x S] cell array of (N x N x T) matrices for each of S subjects                  
%                        Possible values: Unrestricted                                                                      
%                        Default value  : 'MANDATORY INPUT'                                                                 
%                        Input Data Type: any evaluable Matlab expression.                                                  
%                                                                                                                           
% MCMC_InitState:        Object containing initial state of Gibbs sampler                                                   
%                        Input Range  : Unrestricted                                                                        
%                        Default value: n/a                                                                                 
%                        Input Data Type: struct                                                                            
%                                                                                                                           
% NumMcmcIters:          Number of MCMC iterations for spline fitting                                                       
%                        Input Range  : [1  Inf]                                                                            
%                        Default value: 1000                                                                                
%                        Input Data Type: real number (double)                                                              
%                                                                                                                           
% BurnInFractionForMCMC: Fraction of initial MCMC samples to discard (burn in period)                                       
%                        Input Range  : [0        0.99]                                                                     
%                        Default value: 0.5                                                                                 
%                        Input Data Type: real number (double)                                                              
%                                                                                                                           
% ThinningFactor:        Thinning factor for MCMC                                                                           
%                        We keep every kth MCMC sample, where k=ThinningFactor. This is useful when we have limited         
%                        available memory to store MCMC results since successive MCMC estimates are more likely to be       
%                        correlated.                                                                                        
%                        Input Range  : [1  Inf]                                                                            
%                        Default value: 1                                                                                   
%                        Input Data Type: real number (double)                                                              
%                                                                                                                           
% ReturnPosterior:       Return complete posterior distibution in MCMC_LastState                                            
%                        Otherwise return only final state                                                                  
%                        Input Range  : Unrestricted                                                                        
%                        Default value: 1                                                                                   
%                        Input Data Type: boolean                                                                           
%                                                                                                                           
% AppendLastState:       Append new state to initial MCMC state                                                             
%                        Input Range  : Unrestricted                                                                        
%                        Default value: 0                                                                                   
%                        Input Data Type: boolean                                                                           
%                                                                                                                           
% basisCoeffVarPrior:    Variance of basis coefficient gaussian prior                                                       
%                        Larger --> more wiggling allowed                                                                   
%                        Input Range  : [2.220446049250313e-16                    Inf]                                      
%                        Default value: 1000                                                                                
%                        Input Data Type: real number (double)                                                              
%                                                                                                                           
% noiseVarPriorShape:    Shape (D.O.F) of noise variance prior                                                              
%                        This is the "alpha" parameter of the inverse gamma prior distribution. Increasing                  
%                        noiseVarPriorShape --> decreased variance of noise variance distribution.                          
%                        Input Range  : [2.220446049250313e-16                    Inf]                                      
%                        Default value: 0.01                                                                                
%                        Input Data Type: real number (double)                                                              
%                                                                                                                           
% noiseVarPriorScale:    Scale parameter of noise variance prior                                                            
%                        This is the "theta" (1/beta) parameter of inverse gamma prior distribution. Increasing             
%                        noiseVarPriorScale --> right-shift of distribution --> (increase in expected noise variance). In   
%                        general MEAN(noiseVariance) = noiseVarPriorScale/noiseVarPriorShape and MODE(noiseVariance) =      
%                        noiseVarPriorScale/(noiseVarPriorShape-1) for noiseVarPriorShape>=1.                               
%                        Input Range  : [2.220446049250313e-16                    Inf]                                      
%                        Default value: 0.01                                                                                
%                        Input Data Type: real number (double)                                                              
%                                                                                                                           
% TransformData:         Transform data                                                                                     
%                        PreTransFunc is applied before smoothing and its inverse (PostTransFunc) is applied after          
%                        smoothing.                                                                                         
%                        Input Range  : Unrestricted                                                                        
%                        Default value: 1                                                                                   
%                        Input Data Type: boolean                                                                           
%                                                                                                                                                                           
%                                                                                                                           
%     | PreTransFunc:    Pre-smoothing transform function for TSData                                                        
%                        Input Range  : Unrestricted                                                                        
%                        Default value: @(x) log(x+1)                                                                       
%                        Input Data Type: any evaluable Matlab expression.                                                  
%                                                                                                                           
%     | PostTransFunc:   Post-smoothing transform function for TSData                                                       
%                        This should be the inverse of PreTransFunc                                                         
%                        Input Range  : Unrestricted                                                                        
%                        Default value: @(x) exp(x)-1                                                                       
%                        Input Data Type: any evaluable Matlab expression.                                                  
%                                                                                                                           
% VerbosityLevel:        Verbosity level. 0 = no output, 1 = text, 2 = graphical                                            
%                        Possible values: 0,1,2                                                                             
%                        Default value  : 2                                                                                 
%                        Input Data Type: real number (double) 
%
% --------------------------------------------------------------------------------------------------------------------------
% Outputs                                                                                                          
% --------------------------------------------------------------------------------------------------------------------------
%
% FIT                     contains the smooth connectivity estimates
%
% MCMC_LastState          contains the final state of the Gibbs sampler
%                         This is a struct with the fields:
%                             .THETA        [K x Q x NITER]    Q FPCA component vectors, K spline knots
%                             .ALPHA        [Q x R x NITER]    For each of the R connectivity edges in Y, this is the Q-dimensional
%                                                              vector of weights for Q FPCA components. B = mean(ALPHA(q,r,:))
%                                                              is  therefore the expected value of the qth FPCA regression
%                                                              coefficient (Beta-weight) for the rth conn edge)
%                             .ALPHA_BAR    [Q x NITER]        Mean of multivariate gaussian prior for ALPHA (over all edges)
%                             .SIGMA_EPS    [1 x NITER]        Noise variance
%                         where NITER is the number of MCMC iterations (after burn-in)
%                         and additional fields:
%                             .phi_t                      Spline basis functions ([K x T])
%                             .Sigma_theta                THETA cov mat (mvnrnd)
%                             .mu_theta                   THETA mean (mvnrnd)
%                             .Sigma_alpha_i              ALPHA cov mat (mnvnrnd)
%                             .mu_alpha_i                 ALPHA mean (mnvnrnd)
%                             .Sigma_alpha_bar            ALPHA_BAR cov mat (mnvnrnd)
%                             .mu_alpha_bar               ALPHA_BAR mean (mnvnrnd)
%                             .NumSamples                 # MCMC samples stored
% 
%
% Author: Tim Mullen and Wes Thompson, 2010-12, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

g = arg_define([0 Inf],varargin, ...
    arg_norep({'Y','TSData'},mandatory,{},sprintf(['Cell array of data to smooth.\n' ...
              'Generally, Y is a cell array where Y{i} is the T x 1 vector of time-varying (or freq-varying) connectivity for the kth channel pair of the sth subject.\n' ...
              'e.g. Y = {s1(1,1) s1(1,2) s1(1,3) ... s1(N,1) s1(N,2) s1(N,3) ... \n' ...
              '          s2(1,1) s2(1,2) s2(1,3) ... } \n' ...
              'Y can also be a [1 x S] cell array of (N x N x T) matrices for each of S subjects']),'type','expression'), ...
    arg_nogui({'MCMC_InitState'},struct([]),[],'Object containing initial state of Gibbs sampler'), ...
    arg({'nMCMCiters','NumMcmcIters','niters','niter'},1000, [1 Inf], 'Number of MCMC iterations for spline fitting'), ...
    arg({'burnInFraction','BurnInFractionForMCMC'},0.5,[0 0.99],'Fraction of initial MCMC samples to discard (burn in period).'), ...
    arg({'thinFactor','ThinningFactor'},1,[1 Inf],'Thinning factor for MCMC. We keep every kth MCMC sample, where k=ThinningFactor. This is useful when we have limited available memory to store MCMC results since successive MCMC estimates are more likely to be correlated.'), ...
    arg({'returnHistory','ReturnPosterior','ReturnHistory'},true,[],'Return complete posterior distibution in MCMC_LastState. Otherwise return only final state','cat','MCMC State'), ...
    arg({'appendLastState','AppendLastState'},false,[],'Append new state to initial MCMC state','cat','MCMC State'), ... 
    arg({'basisCoeffVarPrior'},1000,[eps Inf],'Variance of basis coefficient gaussian prior. Larger --> more wiggling allowed','cat','Hyperparameters'), ...
    arg({'noiseVarPriorShape'},0.01,[eps Inf],'Shape (D.O.F) of noise variance prior. This is the "alpha" parameter of the inverse gamma prior distribution. Increasing noiseVarPriorShape --> decreased variance of noise variance distribution.','cat','Hyperparameters'), ...
    arg({'noiseVarPriorScale'},0.01,[eps Inf],'Scale parameter of noise variance prior. This is the "theta" (1/beta) parameter of inverse gamma prior distribution. Increasing noiseVarPriorScale --> right-shift of distribution --> (increase in expected noise variance). In general MEAN(noiseVariance) = noiseVarPriorScale/noiseVarPriorShape and MODE(noiseVariance) = noiseVarPriorScale/(noiseVarPriorShape-1) for noiseVarPriorShape>=1.','cat','Hyperparameters'), ...
    arg_subtoggle({'transData','TransformData'},'on',...
    {...
        arg({'preTransFcn','PreTransFunc'},'@(x) log(x+1)',[],'Pre-smoothing transform function for TSData.','type','expression'), ...
        arg({'postTransFcn','PostTransFunc'},'@(x) exp(x)-1',[],'Post-smoothing transform function for TSData. This should be the inverse of PreTransFunc','type','expression') ...
    },'Transform data. PreTransFunc is applied before smoothing and its inverse (PostTransFunc) is applied after smoothing.'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical'));

%     arg_nogui({'logtrans','LogTransform'},true,[],'Log-Transform data before smoothing. Inverse transform is applied after smoothing'), ...

arg_toworkspace(g);

if isempty(MCMC_InitState)
    error('You must provide an initial state for MCMC. See stat_ubsplsm_init');
end

% compute number of burn-in samples
numBurnInSamples = floor(burnInFraction*nMCMCiters);
niterToKeep      = round((nMCMCiters-numBurnInSamples)/thinFactor);
if verb,
    fprintf(hlp_separator());
    fprintf('FPCA Smoothing:\n');
    fprintf('Discarding %d burn-in samples.\n',numBurnInSamples);
    fprintf('Distribution will be thinned by a factor of %d samples.\n',thinFactor);
    fprintf('The distribution of the estimator will have %d samples.\n',niterToKeep);
end

% if Y{i} is in [M x M x T] form, then "flatten" cell entries to [1 x T] form
if ~isvector(Y{1})
    Y = flatten_tensor(Y);
end

% make Y a column vector
Y = Y(:);
% vectorize each element of Y
Y = cellfun(@hlp_vec,Y,'UniformOutput',false);

%% Initialize Parameters
% -------------------------------------------------------------------------
R = size(MCMC_InitState.ALPHA,2);        % number of connectivities to smooth across all subjects
Q = size(MCMC_InitState.ALPHA,1);        % number of FPCA basis functions   
K = size(MCMC_InitState.THETA,1);        % number of smoothing spline knots
T = length(Y{1});                        % dimensionality of connectivity (e.g. number of time points)

% initalize outputs
FIT = zeros(T,R,niterToKeep);            % Smoothed measure estimates
if nargout>1 && returnHistory
    % initialize state object
    MCMC_LastState.THETA        = zeros(K,Q,niterToKeep);      % Q FPCA component vectors, K spline knots
    MCMC_LastState.ALPHA        = zeros(Q,R,niterToKeep);      % For each connectivity edge, the Q-dimensional vector of weights for Q FPCA components
    MCMC_LastState.ALPHA_BAR    = zeros(Q,niterToKeep);        % Mean of multivariate gaussian prior for ALPHA
    MCMC_LastState.SIGMA_EPS    = zeros(1,niterToKeep);        % Noise variance
else
    MCMC_LastState = [];
end
            
% set initial values to last state stored in MCMC object
THETA        = MCMC_InitState.THETA(:,:,end);
ALPHA        = MCMC_InitState.ALPHA(:,:,end);
ALPHA_BAR    = MCMC_InitState.ALPHA_BAR(:,end);
SIGMA_EPS    = MCMC_InitState.SIGMA_EPS(end);
phi_t        = MCMC_InitState.phi_t;  % basis vectors (splines)
if ~appendLastState, clear MCMC_InitState; end

% compute temporal sum of spline basis function covariance matrices:
% e.g. Sigma_phi_sum := sum_t { phi_t(:,t)*phi_t(:,t)' } = phi_t*phi_t';
Sigma_phi_sum = phi_t*phi_t';

%% run MCMC to smooth data
% -------------------------------------------------------------------------
if verb==2
    multiWaitbar('Gibbs sampling','Reset','Color',hlp_getNextUniqueColor);
end

if g.transData.arg_selection
    % pre-transform data
    Y = cellfun(g.transData.preTransFcn,Y,'UniformOutput',false);
end

% constants
% -------------------------------------------------------------------------
Iq              = eye(Q);
% Iq_sp           = speye(Q);
%Iqk             = eye(Q*K);
Sigma_alpha_bar = Iq/(R + 1/basisCoeffVarPrior);
sidx = 0;
for iter=1:nMCMCiters
   

    % Draw ALPHA (spline regression coefficients (weights))
    % ---------------------------------------------------------------------
    
    % compute spline regression coeffs inverse covariance matrix
    Sigma_alpha_i     = Iq + (THETA'*Sigma_phi_sum*THETA)/SIGMA_EPS;
    Sigma_alpha_i     = covfixer(Sigma_alpha_i);
    
    % invert inverse covariance matrix to produce cov mat
    Sigma_alpha_i     = inverse(Sigma_alpha_i);
    Sigma_alpha_i_dbl = double(Sigma_alpha_i);
    % ensure covariance matrix is symmetric and positive definite
    Sigma_alpha_i_dbl = covfixer(Sigma_alpha_i_dbl);
    
    % pre-compute product
    tmpprod = THETA'*phi_t/SIGMA_EPS;                      
    for i=1:R  % for each channel pair
        
        % compute regression coeffs mean
        mu_alpha_i = Sigma_alpha_i*(ALPHA_BAR + tmpprod*Y{i});
        
        % draw spline regression coefficients for this channel pair
        ALPHA(:,i) = mvnrnd(mu_alpha_i,Sigma_alpha_i_dbl)';
    end
    
    
    % Draw ALPHA_BAR (hyperparameter for ALPHA gaussian prior mean)
    % ---------------------------------------------------------------------
    mu_alpha_bar    = sum(ALPHA,2);
    mu_alpha_bar    = Sigma_alpha_bar*mu_alpha_bar;
    ALPHA_BAR       = mvnrnd(mu_alpha_bar,Sigma_alpha_bar)';
    
    % Draw SIGMA_EPS (residual variance a.k.a. noise)
    % ---------------------------------------------------------------------
    pi_eps_shape = noiseVarPriorShape + R*T/2;       % prior distribution shape 
    
    % compute the sum-squared error of residuals (data-fit)^2
    pi_eps_scale = noiseVarPriorScale;          % prior distr. scale param
    tmpprod      = phi_t'*THETA;      % precompute product        
    for i=1:R
        pi_eps_scale = pi_eps_scale ...
               + sum((Y{i}-tmpprod*ALPHA(:,i)).^2 / 2);            
    end
    % draw noise variance estimates from inverse gamma
    % note that we assume the noise variance is identical for all pairs         
    SIGMA_EPS = 1/gamrnd(pi_eps_shape,1/pi_eps_scale);                  
    
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
        Alpha_i    = ALPHA(:,i)';
        Alpha_i_sq = Alpha_i'*Alpha_i;
        for t = 1:T  % for each time window                                     % [!] we should be able to get rid of this loop
            Phi_D_t     = blkdiageye(phi_t(:,t)',Q);                            % [!] this can be pre-computed (for each value of t) and moved to outside the i-loop
            Sigma_theta = Sigma_theta + Phi_D_t'*Alpha_i_sq*Phi_D_t;
            mu_theta    = mu_theta + Phi_D_t'*Alpha_i'*Y{i}(t);
        end
    end
    
    Sigma_theta = Sigma_theta/SIGMA_EPS;
    
    % enforce valid covariance matrix before and after inversion
    Sigma_theta = covfixer(Sigma_theta);
    Sigma_theta = double(inverse(Sigma_theta));
    Sigma_theta = covfixer(Sigma_theta);
    
    % compute mu (mean)
    mu_theta    = mu_theta/SIGMA_EPS;
    mu_theta    = Sigma_theta*mu_theta;
    
    % draw theta
    THETA = reshape(mvnrnd(mu_theta,Sigma_theta)',K,Q);
    
    if verb==2
        multiWaitbar('Gibbs sampling',iter/nMCMCiters);
    end
    
    % Store MCMC estimates of smoothed fits and other parameters
    % ---------------------------------------------------------------------
    if iter>=numBurnInSamples && thinFactor*round(iter/thinFactor)==iter
        sidx = sidx + 1;
        % calculate smoothed fits
        FIT(:,:,sidx) = phi_t'*THETA*ALPHA;
        if nargout > 1 && returnHistory
            % store current MCMC state
            MCMC_LastState.THETA(:,:,sidx)       = THETA;
            MCMC_LastState.ALPHA(:,:,sidx)       = ALPHA;
            MCMC_LastState.ALPHA_BAR(:,sidx)     = ALPHA_BAR;
            MCMC_LastState.SIGMA_EPS(sidx)       = SIGMA_EPS;
        end
    end
end

if g.transData.arg_selection
    % invert transformation of results
    FIT = g.transData.postTransFcn(FIT);
end

if nargout>1
    % finalize MCMC state object
    if ~returnHistory
        % store final values of parameter estimates
        MCMC_LastState.THETA       = THETA;         % Q FPCA component vectors, K spline knots
        MCMC_LastState.ALPHA       = ALPHA;         % Q-dimensional spline regression weights
        MCMC_LastState.ALPHA_BAR   = ALPHA_BAR;     % Gaussian prior mean for ALPHA
        MCMC_LastState.SIGMA_EPS   = SIGMA_EPS;     % Noise variance
    end
    % add additional fields to object
    MCMC_LastState ...
        = catstruct(MCMC_LastState, ...
                    struct( ...
                     'phi_t'            , phi_t,            ...   % Spline basis functions
                     'Sigma_theta'      , Sigma_theta,      ...   % THETA cov mat (mvnrnd)
                     'mu_theta'         , mu_theta,         ...   % THETA mean (mvnrnd)
                     'Sigma_alpha_i'    , Sigma_alpha_i_dbl,...   % ALPHA cov mat (mnvnrnd)
                     'mu_alpha_i'       , mu_alpha_i,       ...   % ALPHA mean (mnvnrnd)
                     'Sigma_alpha_bar'  , Sigma_alpha_bar,  ...   % ALPHA_BAR cov mat (mnvnrnd)
                     'mu_alpha_bar'     , mu_alpha_bar,     ...   % ALPHA_BAR mean (mnvnrnd)
                     'NumSamples'       , fastif(returnHistory,sidx,1), ... % # MCMC samples stored
                     'NumMCMCIters'     , nMCMCiters, ...
                     'initstate'        , false, ...
                     'dataTransform'        , g.transData ...
                     ));
     if appendLastState
        % append current state estimates to initial ones
        MCMC_LastState.THETA       = cat(ndims(THETA)+1,    MCMC_InitState.THETA,     MCMC_LastState.THETA);
        MCMC_LastState.ALPHA       = cat(ndims(ALPHA)+1,    MCMC_InitState.ALPHA,     MCMC_LastState.ALPHA);
        MCMC_LastState.ALPHA_BAR   = cat(ndims(ALPHA_BAR)+1,MCMC_InitState.ALPHA_BAR, MCMC_LastState.ALPHA_BAR);
        MCMC_LastState.SIGMA_EPS   = cat(ndims(SIGMA_EPS)+1,MCMC_InitState.SIGMA_EPS, MCMC_LastState.SIGMA_EPS);
    end
end
     
if verb==2
    % cleanup waitbars
    multiWaitbar('Gibbs sampling'   , 'Close');
    multiWaitbar('Building Splines' , 'Close');
end


% % log transform data
% function [X] = transform(X)
%     X = log(X+1);
% 
% % invert log transform
% function X = untransform(X)
%     X = exp(X)-1;




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
