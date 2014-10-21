% make_eegset()
% --------------------------------------------------------------------------------------------------
function GEEG = grp_make_eegset(varargin)
% build an EEG dataset from provided SIFT data

arg_define([0 Inf],varargin, ...
    arg('EEGref',mandatory,[],'EEG reference dataset. Some fields will be copied from here'), ...
    arg('dipxyz',mandatory,[],'[M x 3] array of XYZ dipole centroids'), ...
    arg('dipcov',[],[],'[3 x 3 x M] dipole covariance matrices'), ...
    arg('dualEquivDipoles',[],[],'List of dipoles that have hemispheric homologues'), ...
    arg('connMat',mandatory,[],'[M x M x T] connectivity matrix'), ...
    arg('connmethod',mandatory,'','Connectivity method name','type','string'), ...
    arg('GEEG',[],[],'Initial group EEG dataset. Conn will be appended to this'), ...
    arg('clusterIds',[],[],'vector of cluster IDs (indices) for each dipole') ...
    );

if isempty(GEEG)
    % build EEG structure
    
    nclusters = size(connMat,1);
    T         = size(connMat,3);
    nbchan    = EEGref.nbchan;
        
    % compute grand means
    % data    = arrayfun(@(EEG_i) mean(EEG_i.data,3)  ,ALLEEG,'UniformOutput',false);
    % icaact  = arrayfun(@(EEG_i) mean(EEG_i.icaact,3),ALLEEG,'UniformOutput',false);
    data = zeros(size(EEGref.data),'single'); %zeros(nbchan,T,'single');
    haveDups  = ~isempty(dualEquivDipoles);
    centroids = [];
    if isempty(clusterIds)
        clusterIds = 1:nclusters; end
    
    % create dipfit structure
    for k=1:nclusters
        % ...store the dipole centroids
        if ismember_bc(k,dualEquivDipoles)
            % (create two dipoles, with homologue on opposite hemisphere)
            centroids(k).posxyz = [dipxyz(k,:) ; dipxyz(k,:).*[1 -1 1]];
            centroids(k).select = [1 2];
        elseif haveDups
            centroids(k).posxyz = [dipxyz(k,:) ; [0 0 0]];
            centroids(k).select = 1;
        else
            centroids(k).posxyz = dipxyz(k,:);
            centroids(k).select = 1;
        end
        % ...store the dipole covariance matrix
        if ~isempty(dipcov)
            centroids(k).covmat = dipcov(:,:,k);
        else
            centroids(k).covmat = [];
        end
        % ...misc fields
        centroids(k).momxyz = eps*ones(2,3);
        centroids(k).rv = 0;
        centroids(k).clustid = clusterIds(k);
    end
    
    % Create final EEG structure
    GEEG                = eeg_emptyset;
    GEEG.data           = data;
    GEEG.icaweights     = eye(nclusters,nbchan);
    GEEG.icasphere      = GEEG.icaweights';
    GEEG.icachansind    = 1:nclusters;
    GEEG.icaact         = [];
    GEEG.srate          = EEGref.srate;
    GEEG.pnts           = EEGref.pnts;
    GEEG.trials         = 1;
    GEEG.chanlocs       = EEGref.chanlocs;
    GEEG.times          = fastif(~isempty(EEGref.times),EEGref.times,((0:(GEEG.pnts-1))/GEEG.srate)*1000);   % ms
    if isempty(GEEG.times), GEEG.times = 0; end
    GEEG.xmin           = GEEG.times(1);
    GEEG.xmax           = GEEG.times(end)/1000;  % sec
    GEEG.nbchan         = EEGref.nbchan;
    GEEG.dipfit         = EEGref.dipfit;
    GEEG.dipfit.model   = centroids;
    GEEG.dipfit.chansel = 1:nclusters;
    
    % check for missing dipfit fields
    if ~isfield(GEEG.dipfit,'mrifile')
        GEEG.dipfit.mrifile = ''; end
    if ~isfield(GEEG.dipfit,'coordformat')
        GEEG.dipfit.coordformat = 'mni'; end
    if ~isfield(GEEG.dipfit,'hdmfile')
        GEEG.dipfit.hdmfile = ''; end
    if ~isfield(GEEG.dipfit,'chanfile')
        GEEG.dipfit.chanfile = []; end 
    
    % intialize SIFT substructures
    GEEG.CAT            = hlp_sift_emptyset(...
                            'curComps',1:nclusters, ...
                            'curComponentNames',strtrim(cellstr(num2str((1:nclusters)')))', ...
                            'nbchan',nclusters, ...
                            'srcdata',GEEG.data);
    GEEG.CAT.PConn      = hlp_sift_emptyconn(...
                            'erWinCenterTimes',EEGref.CAT.Conn.erWinCenterTimes, ...
                            'winCenterTimes',EEGref.CAT.Conn.winCenterTimes, ...
                            'freqs',EEGref.CAT.Conn.freqs, ...
                            'mode','HBMM MCMC');                        
    % check the EEG dataset
    GEEG                = eeg_checkset(GEEG);
end

% Generate final matrix
% insert singleton dimensions if necessary
if length(GEEG.CAT.PConn.winCenterTimes)==1
    connMat = hlp_insertSingletonDim(connMat,4);
end
if length(GEEG.CAT.PConn.freqs)==1
    connMat = hlp_insertSingletonDim(connMat,3);
end

% append this connectivity method
GEEG.CAT.PConn.(connmethod) = connMat;

end