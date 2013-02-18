
% **********************************************************************
% NOTE: THIS FUNCTION IS NOW DEPRECATED AND REPLACED BY stat_surrogate()
% **********************************************************************

function [PConn g] = stat_bootstrap(varargin)

% extract some stuff from inputs for arg defaults
ALLEEG = arg_extract(varargin,'ALLEEG',1);
if ~isempty(ALLEEG)
    if ~isempty(ALLEEG.CAT.Conn)
        Conn = ALLEEG.CAT.Conn(1);
        ConnNames   = hlp_getConnMethodNames(Conn);
        conndef     = ConnNames;
    else
        ConnNames = {''};
        conndef = '';
    end
else
    ConnNames = {''};
    conndef = '';
end
clear ALLEEG Conn;

% get some defaults from ALLEEG
ALLEEG = arg_extract(varargin,'ALLEEG',1);

g = arg_define([0 Inf],varargin, ...
    arg_norep({'ALLEEG','EEG'},mandatory,[],'EEG structure.'), ...
    arg_norep('configs',mandatory,[],'Config superstructure. Contains configuration structures as fields: prepcfg,modfitcfg,conncfg. These are returned from the respective pre_prepData(), est_fitMVAR(), and est_mvarConnectivity() functions'), ...
    arg_subswitch({'mode','Mode'},'Bootstrap', ...
    { ...
    'Bootstrap' { ...
    arg({'nperms','NumPermutations'},200,[],'Number of resamples. This performs efron bootstrapping'), ...
    }, ...
    'Jacknife', { ...
    arg_norep('dummy',[],[]), ...
    }, ...
    'Crossval' { ...
    arg({'nfolds','NumFolds'},10,[],'Number of folds. This performs k-fold resampling'), ...
    } ...
    }, 'Resampling modes. Bootstrap (Efron Bootstrap resampling with replacement), Jacknife (leave-one-out cross-validation), Crossval (k-fold cross-validation)'), ...
    arg_subtoggle({'autosave','AutoSave'},{}, ...
    { ...
    arg({'savefname','FileNamePrefix'},'SIFT_boostrap','','Prefix (including optional path) for autosave file. Data is saved in <FileNamePrefix>.part'), ...
    arg({'savefrq','AutoSaveFrequency'},0.25,[eps 1-eps],'Fractional increment between saves. For example, 0.25 = quarterly saves'), ...
    }, 'Autosave distributions. This will periodically save the computation to a mat file. In case of system crash, bootstrapping can be resumed from this file'), ...
    arg({'connmethods','ConnectivityMethods'},conndef,ConnNames,'Connectivity estimator(s) to bootstrap. All connectivity estimators must be named exactly as they appear in EEG.CAT.Conn','type','logical'), ...
    arg({'verb','VerbosityLevel'},0,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );


arg_toworkspace(g);


% Perform boostrapping
switch lower(mode.arg_selection)
    case 'bootstrap'
        nperms = mode.nperms;
        cv.training = @(x)(randsample(ALLEEG.CAT.trials, ALLEEG.CAT.trials, true));
    case 'jacknife'
        nperms = ALLEEG.CAT.trials;
        cv = cvpartition(ALLEEG.CAT.trials,'leaveout');
    case 'crossval'
        nperms = mode.nfolds;
        cv = cvpartition(ALLEEG.CAT.trials,'kfold',nperms);
end


if verb==1
    progress('init',sprintf('Bootstrapping (%d/%d)...',0,nperms));
    tic
    t0 = toc;
elseif verb==2
    h=statusbar(sprintf('Bootstrapping (x%d)...',nperms));
end

savefrq = round(savefrq*nperms);

for perm=1:nperms
    
    % preprocess data
    EEG = pre_prepData('ALLEEG',ALLEEG,configs.prepcfg,'verb',0,'TrialSubsetToUse',cv.training(perm));
    
    % fit model
    EEG = est_fitMVAR('EEG',EEG,configs.modfitcfg,'verb',0);
    
    % calculate connectivity
    EEG = pop_est_mvarConnectivity(EEG,configs.conncfg,'verb',0);
    
    % append causal estimates to the last dimension of respective arrays in
    % Conn
    if isempty(connmethods)
        connfields = hlp_getConnMethodNames(EEG.CAT.Conn);
    else
        connfields = connmethods;
    end
    
    if perm==1
        for m=1:length(connfields)
            if ~isempty(EEG.CAT.Conn.(connfields{m}))
                PConn.(connfields{m}) = single([]);
            end
        end
    end
    
    for m=1:length(connfields)
        if ~isempty(EEG.CAT.Conn.(connfields{m}))
            dim = ndims(EEG.CAT.Conn.(connfields{m}));
            PConn.(connfields{m}) = cat(dim+1,PConn.(connfields{m}),single(EEG.CAT.Conn.(connfields{m})));
        end
    end
    
    % save checkpoint
    if ~isempty(savefname) && ~mod(perm,savefrq)
        save(sprintf('%s.part',savefname),'PConn','-v7.3');
    end
    
    if verb==1
        te = toc-t0;
        tt = ceil(te*(nperms-perm)/perm);
        progress(perm/nperms,sprintf('Bootstrapping [%d/%d] (EL: %0.1fm / ETA: %0.1fm)',perm,nperms,ceil(te)/60,tt/60));
    elseif verb==2
        statusbar(perm/nperms,h);
    end
    
end

PConn.winCenterTimes = EEG.CAT.Conn.winCenterTimes;
PConn.erWinCenterTimes = EEG.CAT.Conn.erWinCenterTimes;
PConn.freqs = EEG.CAT.Conn.freqs;
PConn.mode = mode.arg_selection;

% clean up
if verb==1
    fprintf('\n')
elseif verb==2
    delete(h);
end

if ~isempty(savefname),
    if exist(sprintf('%s.part',savefname),'file')
        delete(sprintf('%s.part',savefname));
    end
end
