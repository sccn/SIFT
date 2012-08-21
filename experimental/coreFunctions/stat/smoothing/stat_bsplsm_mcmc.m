

function GAMMA_HAT = stat_bsplsm_mcmc(C, bsplmodel,niter,zeta,pi1,pi2)

% C contains [freqs x times] connectivity values for a pair of
% variables/channels
%
% bsplmodel is the b-spline basis model from stat_bsplsm_mkspl
%
% niter is the number of MCMC iterations
%
% niterToKeep is the number of (final) iterations to return in GAMMA_HAT

if nargin<3
    niter = 1000;
end

% commit bspline model to workspace 
arg_toworkspace(bsplmodel);

%% Set up hyperparameters

% prior 
if ~exist('zeta','var')
    zeta=1000; end

% prior
if ~exist('pi1','var')
    pi1=1; end

% prior
if ~exist('pi2','var')
    pi2=1; end


%% Initialize parameters
N=size(Phi_R,1);
NPhiF=size(Phi_F,1);

H=H1*H2;
ETA=zeros(4,niter+1);
DELTA=zeros(H-4,niter+1);
SIGMA_SQ_EPS=zeros(1,niter+1);
SIGMA_SQ_DELTA=zeros(2,niter+1);
GAMMA_HAT=zeros(ntimes,nfreqs,niter+1);

C = C';   % make C be [times x freqs]

% Starting values

ETA(:,1)     = normrnd(0,0.1,[4 1]);
DELTA(:,1)   = normrnd(0,0.1,[H-4 1]);

SIGMA_SQ_EPS(1) = 1;
SIGMA_SQ_DELTA(:,1)  = ones(2,1);


for(iter=1:niter)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               draw ETA                                   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    inv_Sigma_eta=eye(4)/zeta;
    
    mu_eta=sum(reshape(Phi_F,NPhiF,ntimes*nfreqs).*repmat((C(:)'-DELTA(:,iter)'*reshape(Phi_R,N,ntimes*nfreqs))/SIGMA_SQ_EPS(iter),NPhiF,1),2);
    
    for(t=1:ntimes)
        for(f=1:nfreqs)
            inv_Sigma_eta=inv_Sigma_eta+squeeze(Phi_F(:,t,f))*squeeze(Phi_F(:,t,f))'/SIGMA_SQ_EPS(iter);
        end
    end
    
    Sigma_eta=eye(size(inv_Sigma_eta))/inv_Sigma_eta;
    Sigma_eta=.5*(Sigma_eta+Sigma_eta');
    mu_eta=Sigma_eta*mu_eta;
    ETA(:,iter+1)=mvnrnd(mu_eta,Sigma_eta)';
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               draw DELTA                                 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    inv_Sigma_delta=SIGMA_SQ_DELTA(1,iter)\S1_star+SIGMA_SQ_DELTA(2,iter)\S2_star;
    
    mu_delta=sum(reshape(Phi_R,N,ntimes*nfreqs).*repmat((C(:)'-ETA(:,iter+1)'*reshape(Phi_F,NPhiF,ntimes*nfreqs))/SIGMA_SQ_EPS(iter),N,1),2);
    
    for(t=1:ntimes)
        for(f=1:nfreqs)
            inv_Sigma_delta=inv_Sigma_delta+squeeze(Phi_R(:,t,f))*squeeze(Phi_R(:,t,f))'/SIGMA_SQ_EPS(iter);
        end
    end
    Sigma_delta= eye(size(inv_Sigma_delta))/inv_Sigma_delta;
    Sigma_delta=.5*(Sigma_delta+Sigma_delta');
    mu_delta=Sigma_delta*mu_delta;
    DELTA(:,iter+1)=mvnrnd(mu_delta,Sigma_delta)';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               compute GAMMA_HAT                          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    GAMMA_HAT(:,:,iter+1) = reshape(ETA(:,iter+1)'*reshape(Phi_F,NPhiF,ntimes*nfreqs)+DELTA(:,iter+1)'*reshape(Phi_R,N,ntimes*nfreqs),ntimes,nfreqs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               draw SIGMA_SQ_EPS                          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pi1_eps = pi1+.5*nfreqs*ntimes;
    pi2_eps = 0.5*(C-GAMMA_HAT(:,:,iter+1)).^2;
    pi2_eps = sum(pi2_eps(:))+pi2;
    
    
    SIGMA_SQ_EPS(iter+1)=gamrnd(pi1_eps,1/pi2_eps);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               draw SIGMA_SQ_DELTA                        %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pi1_delta=pi1+.5*(H-4);
    pi2_delta=pi2+.5*DELTA(:,iter+1)'*(S1_star+S2_star)*DELTA(:,iter+1);
    SIGMA_SQ_DELTA(:,iter+1)=repmat(gamrnd(pi1_delta,1/pi2_delta),1,2);
    
    
end


% retain only final nitersToKeep of MCMC iterations
nd = ndims(GAMMA_HAT);    
switch nd
    case 3
        GAMMA_HAT = GAMMA_HAT(:,:,end-nitersToKeep+1:end);
    case 2
        GAMMA_HAT = GAMMA_HAT(:,end-nitersToKeep+1:end);
end
            


%%

% smoothY = mean(GAMMA_HAT(:,:,round(iter/2):end),3);
% 
% figure;
% subplot(1,2,1),surf(times,freqs,C'); xlabel('Time'); ylabel('Freq');
% subplot(1,2,2),surf(times,freqs,smoothY'); xlabel('Time'); ylabel('Freq');
% 
% 
% baseline = [1 10];
% basedist = squeeze(mean(GAMMA_HAT(:,baseline,round(iter/2):end),2));
% thresh = prctile(basedist,95,2);
% smoothY(smoothY<repmat(thresh,1,size(smoothY,2))) = 0;
% figure;
% imagesc(times,freqs,smoothY'); set(gca,'Ydir','normal'); xlabel('Time'); ylabel('Freq');
% 


