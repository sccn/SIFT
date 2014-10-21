function [GEEG MCMC_LastState ALLEEG MCMC_InitState] = grp_hbmm(varargin)

% extract some stuff from inputs for arg defaults
ALLEEG = arg_extract(varargin,{'ALLEEG','EEG'},1);

if ~isempty(ALLEEG)
    Conn            = ALLEEG(1).CAT.Conn;
    ConnNames       = hlp_getConnMethodNames(Conn);
else
    Conn = [];
    ConnNames = {''};
end

g = arg_define([0 Inf], varargin, ...
    arg_norep({'ALLEEG','EEG'},mandatory,[],'Array of EEG datasets','type','expression'), ...
    arg_nogui({'MCMC_LastState'},struct([]),[],'Structure containing state of MCMC sampler. This will be used to re-initialize the sampler','type','expression'), ...
    arg({'connmethods','Estimator'},ConnNames{1},ConnNames,'Connectivity estimator to smooth','shape','row'), ...
    arg_subtoggle({'collapseConn','CollapseConn'},{'Conn',Conn,'connmethods',ConnNames,'coldim',{'freq' 'on'}},@hlp_collapseConn,'Select data before smoothing','suppress','connmethods'), ...
    arg_subtoggle({'smoothData','SmoothData'},{},@stat_ubsplsm,'Apply B-spline smoothing to data before group inference','cat','Smoothing','suppress','verb'), ...
    arg({'orthn','Orthnormalize'},true,[],'Orthonormalize smoothing basis functions. This can improve MCMC inference','cat','Smoothing'), ...
    arg_sub({'mcmc_opts','MCMC'},{}, ...
    {...
        arg_sub({'mcmc_init','InitMCMC'},{},@grp_hbmm_init,'Initialize MCMC','suppress',{'verb','transData'}) ...
        arg_sub({'mcmc_run','RunMCMC'},{},@grp_hbmm_mcmc,'MCMC runtime options','suppress',{'verb','transData'}) ...
    },'Perform Markov Chain Monte Carlo inference','cat','MCMC'), ...
    arg_subtoggle({'transData','TransformData'},'on',...
    {...
        arg({'preTransFcn','PreTransFunc'},'@(x) log(x+1)',[],'Pre-MCMC transform function for ConnData.','type','expression'), ...
        arg({'postTransFcn','PostTransFunc'},'@(x) exp(x)-1',[],'Post-MCMC transform function for ConnData. This should be the inverse of PreTransFunc','type','expression') ...
    },'Transform Connectivity data. PreTransFunc is applied before MCMC and its inverse (PostTransFunc) is applied after MCMC. If the data is positive definite, a transform is needed to render it appropriately Gaussian.'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );

% ensure we have a cell array for method names
if ~iscell(g.connmethods)
    g.connmethods = {g.connmethods}; end

% init some stuff
GEEG    = [];
nsubj   = length(ALLEEG);

for si=1:nsubj
    
    Conn = ALLEEG(si).CAT.Conn;
    
    if ~g.collapseConn.arg_selection ...
            && length(Conn.erWinCenterTimes) > 1 ...
            && length(Conn.freqs) > 1
        error('You must collapse at least one dimension (times,freqs). See the "CollapseConn" option.');
    end
    
    % collapse connectivity matrices
    % -------------------------------------------------------------------------
    if g.collapseConn.arg_selection
        Conn = hlp_collapseConn('Conn',Conn,g.collapseConn,'connmethods',g.connmethods);
    end
    
    ALLEEG(si).CAT.Conn = Conn;
end

if length(g.connmethods) > 1
    warning('grp_hbmm:multipleMeasures', ...
        ['You specified more than one connectivity method.' ...
        'This will result in multiple clustering solutions.' ...
        'The output EEG structure accomodates only one dipfit (dipole) structure.' ...
        'We will store only the dipole clustering solution for the first connectivity method.' ...
        'Therefore we recommend instead running this function separately for each connectivity method']);
end

% construct S cell array for dipole locations
% S{i} = [M_i x 3] vector of x,y,z dipole coordinates
% for subject i, for M_i sources
% first identify and strip duplicate dipoles
for si=1:length(ALLEEG)
    dupDipList  = [];
    dipmodel    = ALLEEG(si).dipfit.model;
    for k=1:length(dipmodel)
        % if we have more than one dipole...
        if size(dipmodel(k).posxyz,1) > 1
            % ...check if this is a valid duplicate
            if length(dipmodel(k).select) > 1 || ~all(dipmodel(k).posxyz(2,:)==0)
                % if so, update list of duplicates
                dupDipList = [dupDipList k];
            end
            % ...keep only first dipole
            dipmodel(k).posxyz = dipmodel(k).posxyz(1,:);
        end
    end
    ALLEEG(si).dipfit.model      = dipmodel;
    ALLEEG(si).dipfit.dupDipList = dupDipList;
end

S = arrayfun(@(EEG_i) reshape([EEG_i.dipfit.model.posxyz],3,[])', ...
    ALLEEG,'UniformOutput',false);
    
% smooth desired estimators from Conn structure
% ----------------------------------------------------------------------------------------------
for m=1:length(g.connmethods)
    
    CM = g.connmethods{m};
    
    % determine whether MCMC initial state already exists for this estimator
    NO_MCMC_STATE = isempty(g.MCMC_LastState) || ~isfield(g.MCMC_LastState,CM) || ~isfield(g.MCMC_LastState.(CM),'Z');
    
    % construct B cell array for
    % B{i} = [M_i x M_i x Q] connectivity matrix
    % for subject i, for M_i sources and Q time or frequency points
    B = arrayfun(@(EEG_i) squeeze(EEG_i.CAT.Conn.(CM)), ...
        ALLEEG,'UniformOutput',false);
    
    % pre-transform data
    if g.transData.arg_selection
        B = cellfun(g.transData.preTransFcn,B,'UniformOutput',false);
    end
    
    % smooth the connectivity time-series
    % ----------------------------------------------------------------------------------------------
    if g.smoothData.arg_selection
        
        if isempty(g.MCMC_LastState) || ~isfield(g.MCMC_LastState,CM) || ~isfield(g.MCMC_LastState.(CM), 'UBSPLSM_State')
            % no smoother state detected, we need to smooth the
            % B tensors across last dimension (generally time)...
            
            if g.transData.arg_selection
                % first, disable smoother's pre/post smoothing data transformation
                % (all transformation will be done in the present function)
                g.smoothData.mcmc_opts.mcmc_run.transData.arg_selection = false;
            end
            
            CM = g.connmethods{m};
            switch g.verb
                case 1
                    fprintf(hlp_separator());
                    fprintf('Applying MCMC smoothing to %s...',CM);
                case 2
                    waitbarstr = sprintf('Applying MCMC smoothing to %s',CM);
                    multiWaitbar(waitbarstr,0,'Color',[0.8 0.0 0.1]);
            end
            
            % ...smooth the B tensor estimates
            if nargout > 2
                % return smoothed connectivity distributions
                [FITS, UBSPLSM_State] = stat_ubsplsm(B,g.smoothData);
                if g.transData.arg_selection
                    % invert any transformation applied to B matrix
                    FITS = cellfun(g.transData.postTransFcn,FITS,'UniformOutput',false);
                end
                % store results
                for si=1:length(ALLEEG)
                    ALLEEG(si).CAT.Conn.(CM) = FITS{si};
                end
            else
                [~, UBSPLSM_State] = stat_ubsplsm(B,g.smoothData);
            end
            
            UBSPLSM_DG = UBSPLSM_State.diags;
            UBSPLSM_OD = UBSPLSM_State.offdiags;
            % ...take mean of posterior of smoothing parameters
            if ~isempty(UBSPLSM_DG)
                UBSPLSM_DG.ALPHA     = postmean(UBSPLSM_DG.ALPHA);
                UBSPLSM_DG.THETA     = postmean(UBSPLSM_DG.THETA);
                UBSPLSM_DG.ALPHA_BAR = postmean(UBSPLSM_DG.ALPHA_BAR);
                UBSPLSM_DG.SIGMA_EPS = postmean(UBSPLSM_DG.SIGMA_EPS);
            end
            if ~isempty(UBSPLSM_OD)
                UBSPLSM_OD.ALPHA     = postmean(UBSPLSM_OD.ALPHA);
                UBSPLSM_OD.THETA     = postmean(UBSPLSM_OD.THETA);
                UBSPLSM_OD.ALPHA_BAR = postmean(UBSPLSM_OD.ALPHA_BAR);
                UBSPLSM_OD.SIGMA_EPS = postmean(UBSPLSM_OD.SIGMA_EPS);
            end
            
            % ...orthonormalize mean ALPHA and THETA
            if g.orthn
                if g.verb
                    fprintf('\nOrthonormalizing basis vectors...');
                end
                [dg_inds, od_inds] = deal({});
                for si=1:length(B)
                    % get linear indices of diagonal and off-diagonal elements
                    % of B{si} conn tensor for subject si
                    M_i = size(B{si},1);
                    all_inds = 1:M_i^2;
                    dg_inds{si} = sub2ind([M_i M_i],1:M_i,1:M_i);
                    od_inds{si} = setdiff_bc(all_inds,dg_inds{si});
                end
                % orthonormalize diagonal and off-diagonal coefficients
                UBSPLSM_DG = rot_orthn(UBSPLSM_DG,flatten_tensor(B,dg_inds));
                UBSPLSM_OD = rot_orthn(UBSPLSM_OD,flatten_tensor(B,od_inds));
            end
            
            % ...store the state of the smoother
            UBSPLSM_State.diags     = UBSPLSM_DG;
            UBSPLSM_State.offdiags  = UBSPLSM_OD;
%             MCMC_LastState.(CM).UBSPLSM_State = UBSPLSM_State;
            
            if g.verb
                fprintf('done.\n'); 
            end
            if g.verb==2
                % cleanup waitbar
                multiWaitbar(waitbarstr, 'Close');
            end
        else
            % retrieve the FPCA state
            UBSPLSM_State = g.MCMC_LastState.(CM).UBSPLSM_State;
        end
        
        % replace data B, with smooth basis coefficients (mean posterior of ALPHA)
        % (we will perform MCMC inference on these coefficients)
        smthDG = ~isempty(UBSPLSM_State.diags);
        smthOD = ~isempty(UBSPLSM_State.offdiags);
        
        if smthDG, Q = size(UBSPLSM_State.diags.ALPHA,1);
        else       Q = size(UBSPLSM_State.offdiags.ALPHA,1);
        end
        offset_dg = 0;
        offset_od = 0;
        for si=1:length(B)
            M_i   = size(B{si},1);
            B{si} = zeros(M_i,M_i,Q); % re-init B
            % store diagonals
            if smthDG
                for i=1:M_i
                    B{si}(i,i,:) = UBSPLSM_State.diags.ALPHA(:,i+offset_dg);
                end
                offset_dg = offset_dg + M_i;
            end
            % store off-diagonal
            if smthOD
                k=1;
                for i=1:M_i
                    for j=[1:i-1 i+1:M_i]
                        B{si}(i,j,:) = UBSPLSM_State.offdiags.ALPHA(:,k+offset_od);
                        k=k+1;
                    end
                end
                offset_od = offset_od+M_i^2-M_i;
            end
        end
        
    end % -- smoothData --
    
    
    % Run Heirarchical Bayesian Mixture Model
    % ----------------------------------------------------------------------------------------------
    
    if g.transData.arg_selection
        % first, disable internal pre/post MCMC data transformation
        % (all transformation will be done in the present function)
        g.mcmc_opts.mcmc_run.transData.arg_selection = false;
    end
    
    if NO_MCMC_STATE
        % ...initialize MCMC algorithm
        MCMC_InitState = grp_hbmm_init(B,S,g.mcmc_opts.mcmc_init,'verb',g.verb);
    else
        % ...we have the initial state
        MCMC_InitState = g.MCMC_LastState.(CM);
    end
    
    % ...run MCMC
    MCMC_LastState.(CM) = grp_hbmm_mcmc(B,S,MCMC_InitState,g.mcmc_opts.mcmc_run,'verb',g.verb);
    
    if g.smoothData.arg_selection
        % store the B-spline state again
        MCMC_LastState.(CM).UBSPLSM_State = UBSPLSM_State;
    end

    % Get posterior means of spatial locations & connectivity coefs & indicators
    % ----------------------------------------------------------------------------------------------
    MS = MCMC_LastState.(CM); % abbrev.
    if g.smoothData.arg_selection
        % inference was performed on FPCA beta-weights
        % get group-level smooth connectivity posterior
        % distributions from B-spline basis functions and FPCA weights
        
        have_diags    = ~isempty(MS.UBSPLSM_State.diags);
        have_offdiags = ~isempty(MS.UBSPLSM_State.offdiags);
        
        % get coefficients
        B_BAR = MS.B_BAR;
        
        % get product of [T x Q] spline basis functions
        % and [K x Q] FPCA coefficient matrix (PHI = phi_t'*THETA)
        if have_diags
            PHI_DG =   MS.UBSPLSM_State.diags.phi_t' ...
                     * MS.UBSPLSM_State.diags.THETA;
            T = size(PHI_DG,1); % number of time points
        end
        if have_offdiags
            PHI_OD =   MS.UBSPLSM_State.offdiags.phi_t' ...
                     * MS.UBSPLSM_State.offdiags.THETA;
            T = size(PHI_OD,1); % number of time points
        end
        M = size(B_BAR,1); % number of clusters
        P = size(B_BAR,ndims(B_BAR)); % number of MCMC samples
        
        %       sq = size(PHI); sb = size(B_BAR);
        %       C_BAR = reshape(PHI * B_BAR(:,:),[sq(1) sb(2:end)]);
        
        % initialize data matrix
        C_BAR = zeros(M,M,T,P);
        % reconstruct smoothed connectivities for all MCMC iterations
        for i=1:M
            for j=1:M
                if have_diags && i==j
                    C_BAR(i,j,:,:) = PHI_DG*squish(B_BAR(i,j,:,:));
                elseif have_offdiags && i~=j
                    C_BAR(i,j,:,:) = PHI_OD*squish(B_BAR(i,j,:,:));
                end
            end
        end
    else
        C_BAR = MS.B_BAR;
    end
    
    % invert any transformation applied to B matrix
    if g.transData.arg_selection
        C_BAR = g.transData.postTransFcn(C_BAR);
    end
    
    % get posterior means and covariance of dipole locations
    S_BAR     = postmean(MS.S_BAR);         % mean dipole location
    SIGMA_S   = postmean(MS.SIGMA_S);       % dipole covariance
    %     B_BAR_mean      = postmean(MS.B_BAR);       % parameters mean
    
    % insert connectivity matrix in final EEG dataset
    GEEG = grp_make_eegset('EEGref',ALLEEG(1),'dipxyz',S_BAR,'dipcov',SIGMA_S, ...
                           'connMat',C_BAR,'connmethod',CM,'GEEG',GEEG, ...
                           'dualEquivDipoles',[]);
    % store the number of subjects
    MCMC_LastState.NumSubj = nsubj;
end
end

% ==================================================================================================
% Supplementary Functions
% ==================================================================================================


% rot_orthn()
% --------------------------------------------------------------------------------------------------
function MCMC_State = rot_orthn(MCMC_State,B)
% rotate and scale the basis functions for diagonal connectivities so that the resulting
% smoothing coefficents are orthonormal. This is because the mean
% posterior estimates from smooth_mcmc are not necessarily orthonormal
%
% MCMC_State: contains the final state of the Gibbs sampler
%             This is a struct with the fields:
%             .THETA        [K x Q]    Q FPCA component vectors, K spline knots
%             .ALPHA        [Q x R]    For each of the R connectivity edges in Y, this is the Q-dimensional
%                                      vector of weights for Q FPCA components. B = mean(ALPHA(q,r,:))
%                                      is  therefore the expected value of the qth FPCA regression
%                                      coefficient (Beta-weight) for the rth conn edge)
%             .ALPHA_BAR    [Q x 1]    Mean of multivariate gaussian prior for ALPHA (over all edges)
%             .SIGMA_EPS    [1 x 1]    Noise variance
%
% B  is an [1 x S] cell array of the form
% B = {s1(1,1) s1(1,2) s1(1,3) ... s1(N,1) s1(N,2) s1(N,3) ... s2(1,1) s2(1,2) s2(1,3) ... }
% where sk(i,j) is the Q x 1 vector of connectivity values from channels
% j to i, for the kth subject
B = cellfun(@hlp_vec,B,'UniformOutput',false);
if isempty(MCMC_State), return; end

THETA = MCMC_State.THETA;
ALPHA = MCMC_State.ALPHA;
phi_t = MCMC_State.phi_t;
[K Q] = size(THETA);
R     = size(ALPHA,2);

% compute sample covariance of coefs
S_ALPHA   = R*cov(ALPHA');
% ALPHA_BAR = mean(ALPHA,2);  % mean of posterior estimates of smoothing coefficients
% ALPHA_M   = bsxfun(@minus,ALPHA,ALPHA_BAR);
% S_ALPHA   = ALPHA_M*ALPHA_M';

% orthnonormalize THETA
Sigma_psi    = THETA*S_ALPHA*THETA';
Sigma_psi    = (Sigma_psi+Sigma_psi')/2;
[U D]        = eig(Sigma_psi);
THETA_ORTH   = U(:,(K-Q+1):K)*D((K-Q+1):K,(K-Q+1):K)^.5;
% re-estimate orthonormalized ALPHA coefficients
ALPHA_ORTH   = cellfun(@(B_i) regress(B_i,phi_t'*THETA_ORTH),B,'UniformOutput',false);
ALPHA_ORTH   = cell2mat(ALPHA_ORTH);
S_ALPHA_ORTH = cov(ALPHA_ORTH');

% return orthonormalized theta and alpha
MCMC_State.THETA = THETA_ORTH * sqrtm(S_ALPHA_ORTH);
MCMC_State.ALPHA = S_ALPHA_ORTH^-.5 * ALPHA_ORTH;
% FIT=phi_t'*THETA*ALPHA;

end

% postmean()
% --------------------------------------------------------------------------------------------------
function A = postmean(A)
% return the mean over the last dimension of A
A = mean(A,ndims(A));

end