function [fit_distrib MCMC_LastState]=stat_ubsplsm(Y,K,Q,niters,niterToKeep,verb,smoothDiags,smoothOffDiags,MCMC_LastState)
% univariate b-spline smoothing with fPCA

% Y = cell array of [numChans x numChans x T] connectivities for each subject
% K = number of knots
% Q = number of fpca basis functions

% output:
% fit_distrib: cell array of [numChans x numChans x T x iterations]
%              posterior distribution for each subject

% Author: Tim Mullen and Wes Thompson

if ~exist('verb','var')
    verb = 1;
end

nsubj=length(Y);

% determine the number of channels for each subject
m_i=zeros(nsubj,1);
for i=1:nsubj
    m_i(i)=size(Y{i},1);
end

%% get smoothed time-varying connectivity coefficients

% Diagonals and off-diagonals often have different variances, so we smooth
% them separately

% smooth diagonals
% -------------------------------------------------------------------------
if smoothDiags
    
    if verb==2
        multiWaitbar('Smoothing Diagonals','Reset', ...
                     hlp_getNextUniqueColor('reset'));
    end

    % construct cell array of connectivity matrix
    % diagonal elements for each subject
    CPairs=cell(sum(m_i),1);
    ind=0;
    for i=1:nsubj
        for j=1:m_i(i)
            ind=ind+1;
            CPairs{ind}=squeeze(Y{i}(j,j,:));
        end
    end

    % perform the smoothing
    [fit_diag MCMC_LastState.diag]=stat_ubsplsm_mcmc(CPairs,K,Q,niters,verb);

    % discard burn-in samples
    fit_diag = fit_diag(:,:,end-niterToKeep+1:end);
    
end


% smooth off-diagonals
% -------------------------------------------------------------------------
if smoothOffDiags
    
    if verb
        multiWaitbar('Smoothing Off-Diagonals',1/3);
    end

    % count the total number of pairs
    m_offdiag=sum(m_i.^2-m_i);

    % construct cell array of connectivity matrix
    % off-diagonal elements for each subject
    CPairs=cell(m_offdiag,1);
    ind=0;
    for i=1:nsubj
        for j1=1:m_i(i)
            for j2=1:m_i(i)
                if j1~=j2
                    ind=ind+1;
                    CPairs{ind}=squeeze(Y{i}(j1,j2,:));
                end
            end
        end
    end

    % perform the smoothing
    [alpha alpha_bar theta fit phi_t]=stat_ubsplsm_mcmc(CPairs,K,Q,niters,verb);

    % discard burn-in samples
    fit = fit(:,:,end-niterToKeep+1:end);

end

if verb
    multiwaitbar('Creating final data matrices...',2/3);
end

fit_distrib = cell(1,nsubj);

% store data in [nchs x nchs x time x distribution]
ind=0;
ind_diag = 0;
for i=1:nsubj
    fit_distrib{i} = zeros(m_i(i),m_i(i),size(fit,1),niterToKeep);
    
    for j1=1:m_i(i)
        for j2=1:m_i(i)
            if(not(j1==j2))
                ind=ind+1;
                fit_distrib{i}(j1,j2,:,:)=squeeze(fit(:,ind,:));
            elseif niters~=0
                ind_diag = ind_diag+1;
                fit_distrib{i}(j1,j2,:,:)=squeeze(fit_diag(:,ind_diag,:));
            end
        end
    end
end

if verb
    multiWaitbar('CloseAll');
    pause(0.1);
end

