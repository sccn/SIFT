% hierarchical bayesian mixture model
function [MCMC_State] = grp_hbmm_mcmc(varargin)

% inputs:
%   B: {subj}(M_i x M_i x Q) matrix of (possibly time-varying)  (obtain from EEG(i).CAT.Conn)
%      connectivities or basis coefficients
%   S: {subj}(M_i x 3) matrix of dipole locations (obtain from EEG(i).dipfit)
% automatically determined:
%   N:  number of subjects (length(B))
%   M:  number of clusters (length(N_k))
%   Q:  number of basis coefficients (size(B{1},3))
%
% MCMC_InitState
%     .Z         : {N x 1} CELL array containing [Mi x M] matrices of group indicators for each of N subjects. Z{i}(j,k) = 1 IFF for subject i, source j belongs to cluster k; otherwise 0
%     .S_BAR     : [M x 3] cluster centroid locations
%     .SIGMA_S   : [3 x 3 x M] cluster centroid covariance matrices
%     .B_BAR     : [M x M x Q] group level connectivity means
%     .SIGMA_B   : [M x M] group level connectivity covariance matrix
%     .N_k       : [M x 1] number of subject-level sources that belong to each cluster (group-level source location DOF)
%     .N_k1k2    : [M x M] number of subject-level edges between two clusters (group-level connectivity DOF).

arg_define([0 3], varargin, ...
    arg_norep({'B','Connectivity'},mandatory,[],'Connectivity matrices. For N subjects, this is a cell array B = [1 X N]. B{i} is an (M_i x M_i x Q) matrix of (possibly time-varying) connectivity values or basis coefficients for the ith subject. M_i is the number of components for the ith subject.','type','expression'), ...
    arg_norep({'S','DipoleLocations'},mandatory,[],'Dipole locations. For N subjects, this is a cell array S = [1 X N]. S{i} is an (M_i x 3) matrix of [X Y Z] dipole locations for the ith subject. M_i is the number of components for the ith subject.','type','expression'), ...
    arg_norep({'MCMC_InitState'},mandatory,[],sprintf('Initial state of MCMC sampler. This can be computed from grp_hbmm_initMCMC or can be the output of a previous call to grp_hmbb_mcmc (i.e. the last state of the MCMC iterator). This structure must contain the following fields:\nZ\t: Mi x M matrices of group indicators\nS_BAR\t: cluster centroid locations\nSIGMA_S\t: cluster centroid variances\nB_BAR\t: group level connectivities\nSIGMA_B\t: variances of connectivities\nN_k\t: number of components that belong to each cluster\nN_k1k2\t: pairwise counts of group level clusters')), ...
    arg_sub({'hyperparams','Hyperparameters'},[],...
    {...
        arg({'eta','DipoleLocVarPriorShape'},[],[],'Within-cluster dipole location prior variance D.O.F. Degrees of freedom for inverse-wishart prior distribution for within-cluster dipole location covariance matrices. The a-priori assumption is that dipole location covariance matrices were estimated from eta observations. Leave empty to compute from initial clustering.'), ...
        arg({'SS','DipoleLocVarPriorScale'},[],[],'Within-cluster dipole location prior variance scale. Scale matrix of inverse-wishart prior distribution for within-cluster dipole location covariance matrices. The a-priori assumption is that dipole location covariance matrices were estimated from eta observations with sum of pairwise deviation products SS. Leave empty to compute from initial clustering.'), ...
        arg({'Sm','DipoleLocPriorMean','m'},[],[],'Cluster centroid prior mean hyperparameter. If empty, hyperparameter will be set to initial cluster centroids.'), ...
        arg({'Sr','DipoleLocPriorVar','R'},100*eye(3),[],'Cluster centroid prior covariance hyperparameter. If empty, will use mean covariance of initial clusters.'), ...
        arg({'c','ConnMeanPriorVar'},10000,[eps Inf],'Between-cluster mean connectivity prior variance. Variance of gaussian prior for between-cluster mean connectivity. Larger --> more uncertainty'), ...
        arg({'p1','ConnVarPriorShape','a'},1,[],'Between-cluster mean connectivity variance hyperprior shape. The number of observations (node pairs) is (a priori) assumed to be 2*Shape.'), ...
        arg({'p2','ConnVarPriorScale','b'},1,[],'Between-cluster mean connectivity variance hyperprior scale. The sample variance is (a priori) assumed to be Scale/Shape'), ...
    },'MCMC hyperparameters'), ...
    arg({'nMCMCiters','NumMCMCIters','niter','niters'},1000, [1 Inf], 'Number of MCMC iterations for spline fitting'), ...
    arg({'burnInFraction','BurnInFractionForMCMC'},0.5,[0 0.99],'Fraction of initial MCMC samples to discard (burn in period). The number of MCMC samples is taken to be the smaller of MCMCitersOffDiag and MCMCitersDiag.'), ...
    arg({'thinFactor','ThinningFactor'},1,[1 Inf],'Thinning factor for MCMC. We keep every kth MCMC sample, where k=ThinningFactor. This is useful when we have limited available memory to store MCMC results since successive MCMC estimates are more likely to be correlated.'), ...
    arg_subtoggle({'transData','TransformData'},'on',...
    {...
        arg({'preTransFcn','PreTransFunc'},'@(x) log(x+1)',[],'Pre-MCMC transform function for B.','type','expression'), ...
        arg({'postTransFcn','PostTransFunc'},'@(x) exp(x)-1',[],'Post-MCMC transform function for B. This should be the inverse of PreTransFunc','type','expression') ...
    },'Transform Connectivity data. PreTransFunc is applied before MCMC and its inverse (PostTransFunc) is applied after MCMC. If the data is positive definite, a transform is needed to render it appropriately Gaussian.'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );

%  arg({'normlog','NormalizeAndLogTransform'},true,[],'Transform data before smoothing. Normalize across last dim (time), add 1 and take logarithm. Inverse transform is applied after smoothing'), ...

% set up initial state of MCMC
% --------------------------------------------------------------------------------------------------
if isempty(MCMC_InitState)
    error('You must provide an initial state for the MCMC iterator (MCMC_InitState)');
end

% initialize vars based on last state of MCMC_InitState
% this will initialize the following variables:
% 'Z','S_BAR','SIGMA_S','B_BAR','SIGMA_B','N_k','N_k1k2'
arg_toworkspace(MCMC_InitState);

% define some vars
M   = length(N_k);                      % number of clusters
N   = length(B);                        % number of subjects
M_i = cellfun(@(B_i) size(B_i,1),B);    % number of sources for each subject
Q   = size(B{1},3);                     % connectivity time-series dimension

% initialize prior probabilities of cluster membership
if ~exist('MU','var')
    MU  = ones(M,1)/M; end
% initialize hyperparams
if isempty(hyperparams.eta)
    hyperparams.eta=median(N_k); end
if isempty(hyperparams.SS)
    hyperparams.SS = mean(SIGMA_S,3); end % hyperparams.eta*
% SSinv = double(inverse(hyperparams.SS));
if isempty(hyperparams.Sm)
    hyperparams.Sm = S_BAR; end
if isempty(hyperparams.Sr)
    hyperparams.Sr = hyperparams.SS; end
Sr_inv = double(inverse(hyperparams.Sr)); % SSinv;  %FIXME: setting to SSinv will restore original behavior

% pre-transform data
if transData.arg_selection
    B = cellfun(transData.preTransFunc,B,'UniformOutput',false);
end

% compute number of burn-in samples
numBurnInSamples = floor(burnInFraction*nMCMCiters);
niterToKeep      = round((nMCMCiters-numBurnInSamples)/thinFactor);

if verb,
    fprintf(['I will discard %d burn-in samples.\n' ...
             'I will thin the distribution by a factor of %d samples\n', ...
             'The distribution of the estimator will have %d samples\n'], ...
            numBurnInSamples,thinFactor,niterToKeep);
end

switch verb
    case 1
        fprintf('Running MCMC...');
    case 2
        waitbarstr = sprintf('Running MCMC');
        multiWaitbar(waitbarstr,0,'Color',hlp_getNextUniqueColor());
end

% Run MCMC algorithm
% --------------------------------------------------------------------------------------------------
iter_distrib      = 0;
for iter=1:nMCMCiters
    
    % update waitbar / diagnostics
    switch verb
        case 1
            if ~mod(10,nMCMCiters/iter)
                fprintf('.%0.2f%%',(iter/nMCMCiters)*100);
            end
        case 2
            multiWaitbar(waitbarstr,iter/nMCMCiters);
    end
    
    % Draw S_BAR (group centroid locations)
    % ----------------------------------------------------------------------------------------------
    for k=1:M
        SIGMA_S(:,:,k) = covfixer(SIGMA_S(:,:,k));
        SIGMA_S_inv = double(inverse(squish(SIGMA_S(:,:,k))));
        mu_s_bar_k=zeros(3,1);
        Sigma_s_bar_k = double(inverse(N_k(k)*SIGMA_S_inv+Sr_inv));
        Sigma_s_bar_k = covfixer(Sigma_s_bar_k);
        for i=1:N
            if max(Z{i}(:,k))==1
                j_k=find(Z{i}(:,k)==1);
                for j=1:length(j_k)
                    mu_s_bar_k=mu_s_bar_k+S{i}(j_k(j),:)';
                end
            end
        end
        % OLD VERSION:
%         mu_s_bar_k=Sigma_s_bar_k/SIGMA_S(:,:,k)*mu_s_bar_k;  % FIXME: this is missing the (R^-1)*m term from step 1 of alg.
        mu_s_bar_k=Sigma_s_bar_k*(Sr_inv*hyperparams.Sm(k,:)' + SIGMA_S_inv*mu_s_bar_k);
        S_BAR(k,:)=mvnrnd(mu_s_bar_k,Sigma_s_bar_k);
    end

    
    % Draw SIGMA_S (group centroid covariance matrices)
    % ----------------------------------------------------------------------------------------------
    for k=1:M
        DIV = 0.5;
        eta_k=hyperparams.eta+DIV*N_k(k);
        SS_k=hyperparams.SS;  % FIXME: this is the original version
        for i=1:N
            if max(Z{i}(:,k))==1
                j_k=find(Z{i}(:,k)==1);
                for j=1:length(j_k)
                    SS_k=SS_k+DIV*(S{i}(j_k(j),:)-S_BAR(k,:))'*(S{i}(j_k(j),:)-S_BAR(k,:));
                end
            end
        end
        inv_SS_k      = covfixer(double(inverse(covfixer(SS_k))));
        SIGMA_S(:,:,k)= double(inverse(wishrnd(inv_SS_k,eta_k)));
    end  
    
    % Draw B_BAR (group mean connectivites)
    % ----------------------------------------------------------------------------------------------
    for k1=1:M
        for k2=1:M
            mu_b_bar_k1k2=zeros(Q,1);
            Sigma_sq_b_k1k2=1/(1/hyperparams.c+N_k1k2(k1,k2)/SIGMA_B(k1,k2))*eye(Q);
            for i=1:N
                n_k1=sum(Z{i}(:,k1));
                n_k2=sum(Z{i}(:,k2));
                if(and(n_k1>=1,n_k2>=1))
                    j_k1=find(Z{i}(:,k1)==1);
                    j_k2=find(Z{i}(:,k2)==1);
                    for j1=1:n_k1
                        for j2=1:n_k2
                            mu_b_bar_k1k2=mu_b_bar_k1k2+squish(B{i}(j_k1(j1),j_k2(j2),:));
                        end
                    end
                end
            end
            mu_b_bar_k1k2=Sigma_sq_b_k1k2*mu_b_bar_k1k2/SIGMA_B(k1,k2);
            B_BAR(k1,k2,:)=mvnrnd(mu_b_bar_k1k2,Sigma_sq_b_k1k2);
        end
    end
    
    % Draw SIGMA_B (group connectivity variances)
    % ----------------------------------------------------------------------------------------------
    for k1=1:M
        for k2=1:M
            DIV = quickif(k1==k2,0.5,2); % auto-connections have 1/2 the d.o.f.
            p1_k1k2=hyperparams.p1+DIV*Q*N_k1k2(k1,k2);
            p2_k1k2=hyperparams.p2;
            for i=1:N
                n_k1=sum(Z{i}(:,k1));
                n_k2=sum(Z{i}(:,k2));
                if(and(n_k1>=1,n_k2>=1))
                    j_k1=find(Z{i}(:,k1)==1);
                    j_k2=find(Z{i}(:,k2)==1);
                    for j1=1:n_k1
                        for j2=1:n_k2
                            p2_k1k2=p2_k1k2+DIV*sum((B{i}(j_k1(j1),j_k2(j2),:)-B_BAR(k1,k2,:)).^2);
                        end
                    end
                end
            end
            SIGMA_B(k1,k2)=1/gamrnd(p1_k1k2,1/p2_k1k2);
        end
    end


    % Draw Z (indicators of group membership) 
    % ----------------------------------------------------------------------------------------------
    ind=randperm(N);
    for i=ind
        ind_i=randperm(M_i(i));
        for j=ind_i
            Log_p_ij=zeros(M,1);
            for k=1:M
                Sigma_s_ijk=SIGMA_S(:,:,k);
                s_ijk=S{i}(j,:)-S_BAR(k,:);
                log_p_s_ijk=-.5*log(det(Sigma_s_ijk))-.5*s_ijk/Sigma_s_ijk*s_ijk';
                log_p_b_ijk=0;
                for j1=1:M_i(i)
                    if j==j1
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(k,k)))...
                        -.5*sum((B{i}(j,j,:)-B_BAR(k,k,:)).^2)/SIGMA_B(k,k);
                    end
                    if and(not(j==j1),find(Z{i}(j1,:)==1)==k)
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(k,k)))...
                        -.5*sum((B{i}(j1,j1,:)-B_BAR(Z{i}(j1,:)==1,Z{i}(j1,:)==1,:)).^2)/...
                            SIGMA_B(Z{i}(j,:)==1,Z{i}(j,:)==1);                    
                    end
                    if not(find(Z{i}(j,:)==1)==find(Z{i}(j1,:)==1))
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(k,Z{i}(j1,:)==1)))...
                        -.5*sum((B{i}(j,j1,:)-B_BAR(k,Z{i}(j1,:)==1,:)).^2)/...
                              SIGMA_B(k,Z{i}(j1,:)==1);
                    end
                end
                for j2=1:M_i(i)
                    if not(find(Z{i}(j2,:)==1)==find(Z{i}(j,:)==1))
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(Z{i}(j2,:)==1,k)))...
                        -.5*sum((B{i}(j2,j,:)-B_BAR(Z{i}(j2,:)==1,k,:)).^2)/...
                                SIGMA_B(Z{i}(j2,:)==1,k);
                    end
                end
                Log_p_ij(k)=MU(k)+log_p_s_ijk+log_p_b_ijk;
            end
            Log_p_ij=Log_p_ij-max(Log_p_ij);
            P_ij=exp(Log_p_ij)/sum(exp(Log_p_ij));
            v=rand;
            Z{i}(j,:)=zeros(1,size(Z{i},2)); %0*Z{i}(j,:);
            if v<=P_ij(1)
                   Z{i}(j,1)=1;
                   log_p_ij=Log_p_ij(1);
            end
            for k=1:(M-1)
               if(and(v>sum(P_ij(1:k)),v<=sum(P_ij(1:(k+1)))))
                   Z{i}(j,k+1)=1;
                   log_p_ij=Log_p_ij(k);
               end
            end
        end
    end

    % Draw MU (probabilities of clusters)   
    % ----------------------------------------------------------------------------------------------
    N_k1k2=zeros(M);
    for k1=1:M
        for k2=1:M
            for i=1:N
                n_k1=sum(Z{i}(:,k1));
                n_k2=sum(Z{i}(:,k2));
                if(and(n_k1>=1,n_k2>=1))
                    if not(k1==k2)
                        N_k1k2(k1,k2)=N_k1k2(k1,k2)+n_k1*n_k2;
                    end                
                    if k1==k2
                        N_k1k2(k1,k2)=N_k1k2(k1,k2)+n_k1;
                    end
                end
            end
        end
    end
    N_k=diag(N_k1k2);
    MU=drchrnd(1/M*ones(1,M)+N_k',1);

    
    % Store posterior draws  
    % ----------------------------------------------------------------------------------------------
    if iter>=numBurnInSamples && thinFactor*round(iter/thinFactor)==iter
        iter_distrib=iter_distrib+1;
        
        if iter_distrib == 1
            % initialize arrays of posterior draws
            Z_distrib         = cell(N, 1);
            for si=1:N
                Z_distrib{si} = false(M_i(si),M,niterToKeep); 
            end
            S_BAR_distrib     = zeros([niterToKeep,   size(S_BAR)  ], class(S_BAR));
            SIGMA_S_distrib   = zeros([niterToKeep,   size(SIGMA_S)], class(SIGMA_S));
            B_BAR_distrib     = zeros([niterToKeep,   size(B_BAR)  ], class(B_BAR));
            SIGMA_B_distrib   = zeros([niterToKeep,   size(SIGMA_B)], class(SIGMA_B));
            N_k_distrib       = zeros([niterToKeep,   size(N_k)    ], class(N_k));
            N_k1k2_distrib    = zeros([niterToKeep,   size(N_k1k2) ], class(N_k1k2));
            MU_distrib        = zeros([niterToKeep,   size(MU)     ], class(MU));
        end
        
        % store results for this iteration 
        % (extra trailing colons are to handle any number of dimensions)
        for si=1:N
        	Z_distrib{si}(:,:,iter_distrib) = Z{si}; 
        end
        S_BAR_distrib(iter_distrib,:,:,:,:)     = S_BAR;
        SIGMA_S_distrib(iter_distrib,:,:,:,:)   = SIGMA_S; % ./repmat(reshape(N_k,[1 1 length(N_k)]), [size(SIGMA_S,1) size(SIGMA_S,2)])
        B_BAR_distrib(iter_distrib,:,:,:,:)     = B_BAR;
        SIGMA_B_distrib(iter_distrib,:,:,:,:)   = SIGMA_B;
        N_k_distrib(iter_distrib,:,:,:,:)       = N_k;
        N_k1k2_distrib(iter_distrib,:,:,:,:)    = N_k1k2;
        MU_distrib(iter_distrib,:,:,:,:)        = MU;
    end
end

if transData.arg_selection
    % invert transformation, applied to B matrix
    B_BAR_distrib = transData.postTransFcn(B_BAR_distrib);
end

switch verb
    case 1
        fprintf('done.\n');
    case 2
        % cleanup waitbar
        multiWaitbar(waitbarstr, 'Close');
end
            
% construct MCMC state object
% --------------------------------------------------------------------------------------------------
% MCMC samples (draws) are stored in last dimenions of each array
varnames = {'Z_distrib',        ...
            'S_BAR_distrib',    ...
            'SIGMA_S_distrib',  ...
            'B_BAR_distrib',    ...
            'SIGMA_B_distrib',  ...
            'MU_distrib',       ...
            'N_k1k2_distrib',   ...
            'N_k_distrib'};
for i = 1:length(varnames)
    % strip trailing '_distrib' from state fieldnames
    fname = strrep(varnames{i},'_distrib','');
    % store state estimate with MCMC samples in last dimension
    A  = eval(varnames{i}); 
    nd = ndims(A);
    MCMC_State.(fname) = permute(A,[2:nd 1]);
end
MCMC_State.initstate = false;


% ==================================================================================================
% Helper functions
% ==================================================================================================


function r = drchrnd(a,n)
% sample from a dirichlet distribution
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);
