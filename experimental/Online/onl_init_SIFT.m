function EEG = onl_init_SIFT(EEG,curComps,curComponentNames,usechannels)

    if nargin<4
        usechannels = false;
    end
    
    EEG.CAT.MODEL               = [];
    EEG.CAT.curComps            = curComps;
    EEG.CAT.curComponentNames   = curComponentNames;
    EEG.CAT.nbchan              = length(curComps);
    EEG.times                   = linspace(EEG.xmin*1000,EEG.xmax*1000,EEG.pnts);
    
    if usechannels
        
        % use identity matrices for ICA solution
        [EEG.icaweights EEG.icasphere EEG.icawinv] = deal(eye(size(EEG.data,1)));
        
        % copy data into EEG.icaact and EEG.CAT.srcdata
        [EEG.CAT.srcdata EEG.icaact] = deal(EEG.data(curComps,:));

        % create dipfit structure for plotting 3D channel locations
        % (NOTE: this should be done in run_readmindo)
        EEG.dipfit.hdmfile  = '';
        EEG.dipfit.mrifile  = '';
        EEG.dipfit.chanfile = '';
        EEG.dipfit.chansel  = 1:EEG.nbchan;
        EEG.dipfit.coordformat = 'Spherical';
        EEG.dipfit.coord_transform = [0 0 0 0 0 0 804.9000 1025 784.9000];
                
        for i=1:EEG.nbchan
            EEG.dipfit.model(i).posxyz   = [EEG.chanlocs(i).X EEG.chanlocs(i).Y EEG.chanlocs(i).Z]*900;
            EEG.dipfit.model(i).momxyz  = [0.1 0.1 0.1];
            EEG.dipfit.model(i).rv      = 0.1;
            EEG.dipfit.model(i).select  = 1;
            EEG.dipfit.model(i).diffmap = [];
            EEG.dipfit.model(i).sourcepot = [];
            EEG.dipfit.model(i).datapot  = [];
        end
    else
        
        [EEG.CAT.srcdata EEG.icaact] = deal(EEG.icaweights(curComps,:)*EEG.data);
        
    end
   