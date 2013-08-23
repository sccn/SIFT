% TestDataPath = '/Users/timmullen/Documents/WORK/SIFT/Code/SIFT_bitbucket/experimental/Online/Mindo/mindo_test.set/mindo_test.set';
% TestDataPath = 'C:\Users\tim\Documents\MainlyMozart\SIFT_bitbucket\experimental\Online\eb72_continuous.set';
% TestDataPath = '/Users/timmullen/Documents/WORK/SIFT/Code/SIFT_bitbucket/experimental/Online/MikeChi/MikeChi_TestData.set';

% USE THIS FOR LIVE EVENT:
% this is the path to the 'calibration set'
TestDataPath = 'C:\Users\christian\Documents\MainlyMozart\SIFT_bitbucket\experimental\Online\DataSets\newrecording3.set';  
 % this is the path to the dataset with the pre-computed ica weight matrices, dipfit, etc  
ICAOfflineDataPath = 'C:\Users\christian\Documents\MainlyMozart\SIFT_bitbucket\experimental\Online\DataSets\338_020212_Speech_tweaked.set';  %'C:\Users\tim\Documents\MainlyMozart\SIFT_bitbucket\experimental\Online\DATASET_WITH_ICA_TESTING.set'; 
% This is the path to brainmovie configuration file with variable 'BMCFG'
BMCONFIG_FILE = 'C:\Users\christian\Documents\MainlyMozart\SIFT_bitbucket\experimental\Online\MM_BMCFG.mat'; %'/Users/timmullen/Documents/WORK/SIFT/Code/SIFT_bitbucket/experimental/Online/BMConfigs/MINDO16_BMCFG.mat';

%% load dataset and get ICA and dipole stuff
tmp = pop_loadset(ICAOfflineDataPath);
OfflineDataset.icawinv          = tmp.icawinv;
OfflineDataset.icaweights       = tmp.icaweights;
OfflineDataset.icasphere        = tmp.icasphere;
OfflineDataset.icachansind      = tmp.icachansind;
OfflineDataset.dipfit           = tmp.dipfit;

% get fieldnames
OfflineDatasetFieldnames = fieldnames(OfflineDataset);

% clear the tmp dataset
clear tmp;

%% PARAMETERS
UsingChristianBIOSEMI = true;                   % set this when we are loading datasets recorded from christian's biosemi interface

Components_Speech_Super = [4 7 11 16 23 25 31 35 38]; % 55 62 9 17 19 14];
Components_Speech_Good = [10 21];
Components_Speech_Questionable = [8 26 41 48 56 64];

% -------- PARAMS -----------
Channels = [Components_Speech_Super]; % Components_Speech_Good Components_Speech_Questionable];                  % these are the components/channels to which we'll fit our multivariate model  % eb72: [8 11 12 13 14 15]
StreamViewerChannels = Channels;
WindowLengthSec = 5;                            % sliding window length in seconds
WindowStepSizeSec = WindowLengthSec;            % sliding window step size in seconds (ignore for online)
Frequencies = 1:30;
ModelOrder = 10;
PowerChannels = Channels;                       % channel for which to plot power
ConnectivityMeasures = {'nPDC','S'};            % see doc est_mvtransfer for a full list of codes
SamplingRate = 128;

BENCHMARK = false;
CHECK_WHITENESS = false;
PLOTDATA = true;
PLOTSPECTRA = false;
USE_CHANNELS = false;
REREF = false;
% ReferenceChannels = [];
% ReferenceExclusions = [147:152];

% Set up bad channels (detect this later...)
ChansToSelect = 1:128;
BadChannels = zeros(1,length(ChansToSelect)); % init
ForceBadChannels = []; % Indices of channels to force to bad
DataRejectStdevFactorThresholds = [1/5 5];

% parameters for Edge/Node Size/Color limits adaptation
ADAPT_LIMITS = true;
ADAPTATION_HL = 10; % half-life of moving average (in frames)
MEMFACTOR = 2/((ADAPTATION_HL * 2.8854)+1);
UPDATE_INTERVAL = 1;

FilterParams = [1 2 40 50]; %[0.5 1]; %[0.5 1 6 10];
FilterType = 'bandpass';
% FilterType = 'lowpass';

%% Set up the automatic artifact rejection system
if ~isempty(DataRejectStdevFactorThresholds)
    fprintf('Computing automatic artifact rejection thresholds...\n');
    tmpEEG = pop_loadset(ICAOfflineDataPath);
    % get the upper and lower percentiles of the variance distribution of the data, which we'll use for
    % channel rejection.
    stdevs = std(tmpEEG.data(ChansToSelect,:),0,2);  %
    DataRejectionThresholds = [DataRejectStdevFactorThresholds(1)*min(stdevs) DataRejectStdevFactorThresholds(2)*max(stdevs)];
%     DataRejectionThresholds = prctile(var(tmp.EEG.data(ChansToSelect,:),0,2),[DataRejectStdevFactorThresholds(1) DataRejectStdevFactorThresholds(2)]); 
    clear tmpEEG;
end

%% SET UP FILTERS
% load calibration set
EEG = exp_eval(io_loadset(ICAOfflineDataPath));   % TestDataPath
processed = EEG;
% re-reference
% if REREF
%     processed = exp_eval(flt_reref('Signal',processed,'ReferenceChannels',ReferenceChannels,'ExcludeChannels',ReferenceExclusions,'KeepReference',true));
% end

if ~isempty(ChansToSelect)
    labels = {EEG.chanlocs.labels};
    processed = exp_eval(flt_selchans(processed,labels(ChansToSelect)));
end

% band-pass filter
if ~isempty(FilterParams)
    processed = exp_eval(flt_iir(processed,FilterParams,FilterType));
end


%% RUN PRE-RECORDED DATA
rawdata = io_loadset(ICAOfflineDataPath);
run_readdataset('mystream',rawdata,25);

%% run biosemi system:
% run_readlsl('MatlabStream','mystream','SelectionValue','EEG');

%% Run this to write the current stream to file
% run_writedataset('SourceStream','mystream');

%% put a pipeline on top of the data stream that replicates the processing applied to processed and continues it on new data
pipln = onl_newpipeline(processed,{'mystream'});

%% clear temporary data
% clear processed EEG

%% -------------------------------------

pause(1);

ChannelIDs = strtrim(cellstr(num2str(Channels')))';  % convert list of components to cell array of strings
% ChannelIDs = {EEG.chanlocs(Channels).labels};
PowerChannelsIdx = find(ismember_bc(Channels,PowerChannels));

% reset the state
clear mvar_glADMM;
%
try 
    tmp = load(BMCONFIG_FILE);
    BMCFG = tmp.BMCFG;
    BMCFG.connmethod = ConnectivityMeasures{1};
    BMCFG.BMopts.bmopts_suppl = {'title',{'Multivariate ','Granger Causality  '}};
    BMCFG.BMopts.caption = true;
catch
    BMCFG = struct([]);
end
closed = false;
init = true;
mode = 'init_and_render';


if PLOTDATA
    % note: replace 'basewsDatasetName','' to view all channels
    streamViewerHandle = vis_dataStreamViewer('streamname','mystream','basewsDatasetName','EEG','channelsToDisplay',StreamViewerChannels,'spacing',1, ...
        'closeRequestFcn','svtimer = timerfind(''Tag'',''StreamViewerTimer''); stop(svtimer); delete(svtimer); delete(gcbf);');
    dataTimer = timer('Tag','StreamViewerTimer','StartDelay',0.2, 'Period', 1, 'ExecutionMode','fixedRate', ...
        'TimerFcn',@(x,y) vis_dataStreamViewer('streamname','mystream','basewsDatasetName','EEG','axisHandle',streamViewerHandle,'draw',true,'channelsToDisplay',StreamViewerChannels));  %'StreamViewerTimer'
    start(dataTimer);
end

if PLOTSPECTRA
    hpower = figure('DeleteFcn',@(varargin)evalin('base','closed=true;'));
    haxpwr = axes;
    pwr = zeros(length(PowerChannelsIdx),length(Frequencies));
end

counter = 0;

numberOfLimitsCalculationsSoFar = 0;
numberOfRunsSofar = 0;

BadChannels(ForceBadChannels) = 1;

% poll the buffer periodically
while ~closed
    
    % grab the last WindowLengthSec chunk of data from the stream
    [EEG,pipln] = onl_filtered(pipln, EEG.srate*WindowLengthSec);
%     EEG = onl_peek('mystream',WindowLengthSec);

    if EEG.pnts < EEG.srate*WindowLengthSec
        disp('No Data!');
        continue;
    end

    if ~isempty(DataRejectStdevFactorThresholds)
        % reject channels with variance outside desired bounds
        stdevs = std(EEG.data,0,2);
        BadChannels = (stdevs < DataRejectionThresholds(1) | stdevs > DataRejectionThresholds(2));
        
        BadChannels(ForceBadChannels) = 1;
        
        if all(BadChannels)
            disp('All channels are bad!');
            continue;
        end
    end
    goodchannels = ~BadChannels;
        
    if REREF
        % do common average re-referencing
        EEG = flt_comAvgReref(EEG,goodchannels);
    end
    
    for fn = 1:length(OfflineDatasetFieldnames)
        % insert the stuff from Offline dataset (icawinv, dipfit, etc)
        EEG.(OfflineDatasetFieldnames{fn}) = OfflineDataset.(OfflineDatasetFieldnames{fn});
    end
    
    % initialize the SIFT datastructures
    EEG = onl_init_SIFT(EEG,Channels,ChannelIDs,USE_CHANNELS,goodchannels);
    
    % just a minor fix
    if UsingChristianBIOSEMI
        EEG.chaninfo.nosedir = '+Y';
    end
    
    tic
    
    % fit MVAR model
    try
        EEG = pop_est_fitMVAR(EEG,0,'normalize',{'time'},'detrend','linear','algorithm','gladmm','morder',ModelOrder, ...
            'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'epochTimeLims',[EEG.xmin EEG.xmax],'verb',0, ...
            'gladmm',struct(...
            'lambda',0.0005,        ...   %0.0013
            'alpha',1.8,            ...
            'rho',2,                ...
            'verb',0,               ...
            'max_iter',500,         ...
            'rho_update',false,     ...
            'RhoUpdateIncr',1.1,    ...
            'RhoUpdateDecr',1.1,    ...
            'lambda_update',true,   ...
            'lambda_update_thresh',10^-5,   ...
            'lambda_update_factor',1.1,     ...
            'lambda_update_count',5)        ...
            );
    catch
        disp('Bad Data: ADMM did not converge!');
        continue;
    end
    
    % get connectivity
    EEG.CAT.Conn = est_mvarConnectivity(EEG,EEG.CAT.MODEL,'connmethods',ConnectivityMeasures,'freqs',Frequencies,'verb',0);
    
    EEG.CAT.Conn = hlp_absvalsq(EEG.CAT.Conn,ConnectivityMeasures,false,false);
    
    if isfield(EEG.CAT.Conn,'S')
%         EEG.CAT.Conn = hlp_absvalsq(EEG.CAT.Conn,{'S'});
        EEG.CAT.Conn.S = 10*log10(EEG.CAT.Conn.S);
    end
    
    tfit = toc;
    
    if init
        % initialize brainmovie
        if ~isempty(BMCFG)
            BMCFG.showNodeLabels.nodelabels = ChannelIDs;
            hbmcfg = gui_causalBrainMovie3D_online(EEG,EEG.CAT.Conn,struct('arg_direct',0),BMCFG);
        else
            hbmcfg = gui_causalBrainMovie3D_online(EEG,EEG.CAT.Conn,struct('arg_direct',0),'showNodeLabels',{'nodelabels',ChannelIDs});
        end
        waitfor(hbmcfg,'UserData','init');
        %         BMCFG.BMopts.speedy = true;
        bmoptschanged = true;
        init = false;
        figh = [];
        bmvars = [];
    end
    
    if bmoptschanged
        mode = 'init_and_render';
        bmoptschanged = false;
        bmvars = [];
        
        % reset the limits
        numberOfRunsSofar = 0;
        numberOfLimitsCalculationsSoFar = 0;
    end
    
    tic
    
    if numberOfRunsSofar > 0
        % update the data limits
        BMCFG.BMopts.graphColorAndScaling.nodeSizeDataRange  = [lastNodeSizeMin lastNodeSizeMax];
        BMCFG.BMopts.graphColorAndScaling.edgeSizeDataRange  = [lastEdgeSizeMin lastEdgeSizeMax];
        BMCFG.BMopts.graphColorAndScaling.nodeColorDataRange = [lastNodeColorMin lastNodeColorMax];
        BMCFG.BMopts.graphColorAndScaling.edgeColorDataRange = [lastEdgeColorMin lastEdgeColorMax];
    end
    
    opts = hlp_mergeVarargin(BMCFG.BMopts,'figurehandle',figh,'mode',mode,'speedy',true,'vars',bmvars);
        
    [tmp handles BMout] = vis_causalBrainMovie3D(EEG,EEG.CAT.Conn,BMCFG,'timeRange',[],'BMopts',opts);
    
    tdraw = toc;
    
    if ADAPT_LIMITS
        
        if (numberOfRunsSofar < 10 || mod(numberOfRunsSofar,UPDATE_INTERVAL) == 0) && numberOfRunsSofar < Inf % 50
            
            if numberOfRunsSofar < 10
                
                lastNodeSizeMin  =  min(BMout.NodeSize);
                lastNodeSizeMax  =  max(BMout.NodeSize);
                
                lastNodeColorMin =  min(BMout.NodeColor);
                lastNodeColorMax =  max(BMout.NodeColor);
                
                lastEdgeSizeMin  =  min(BMout.EdgeSize(:));
                lastEdgeSizeMax  =  max(BMout.EdgeSize(:));
                
                lastEdgeColorMin =  min(BMout.EdgeColor(:));
                lastEdgeColorMax =  max(BMout.EdgeColor(:));
                
            else
                lastNodeSizeMin  = MEMFACTOR * min(BMout.NodeSize) + (1-MEMFACTOR) * lastNodeSizeMin; % (MEMFACTOR*(lastNodeSizeMin * numberOfLimitsCalculationsSoFar) + min(BMout.NodeSize)) / (numberOfLimitsCalculationsSoFar+1);
                lastNodeSizeMax  = MEMFACTOR * max(BMout.NodeSize) + (1-MEMFACTOR) * lastNodeSizeMax; %(MEMFACTOR*(lastNodeSizeMax * numberOfLimitsCalculationsSoFar) + max(BMout.NodeSize)) / (numberOfLimitsCalculationsSoFar+1);
                
                lastNodeColorMin = MEMFACTOR * min(BMout.NodeColor) + (1-MEMFACTOR) * lastNodeColorMin; %(MEMFACTOR*(lastNodeColorMin * numberOfLimitsCalculationsSoFar) + min(BMout.NodeColor)) / (numberOfLimitsCalculationsSoFar+1);
                lastNodeColorMax = MEMFACTOR * max(BMout.NodeColor) + (1-MEMFACTOR) * lastNodeColorMax; %(MEMFACTOR*(lastNodeColorMax * numberOfLimitsCalculationsSoFar) + max(BMout.NodeColor)) / (numberOfLimitsCalculationsSoFar+1);
                
                lastEdgeSizeMin  = MEMFACTOR * min(BMout.EdgeSize(:)) + (1-MEMFACTOR) * lastEdgeSizeMin; %(MEMFACTOR*(lastEdgeSizeMin * numberOfLimitsCalculationsSoFar) + min(BMout.EdgeSize(:))) / (numberOfLimitsCalculationsSoFar+1);
                lastEdgeSizeMax  = MEMFACTOR * max(BMout.EdgeSize(:)) + (1-MEMFACTOR) * lastEdgeSizeMax; %(MEMFACTOR*(lastEdgeSizeMax * numberOfLimitsCalculationsSoFar) + max(BMout.EdgeSize(:))) / (numberOfLimitsCalculationsSoFar+1);
                
                lastEdgeColorMin = MEMFACTOR * min(BMout.EdgeSize(:)) + (1-MEMFACTOR) * lastEdgeColorMin;  %(MEMFACTOR*(lastEdgeColorMin * numberOfLimitsCalculationsSoFar) + min(BMout.EdgeColor(:))) / (numberOfLimitsCalculationsSoFar+1);
                lastEdgeColorMax = MEMFACTOR * max(BMout.EdgeSize(:)) + (1-MEMFACTOR) * lastEdgeColorMax;  %(MEMFACTOR*(lastEdgeColorMax * numberOfLimitsCalculationsSoFar) + max(BMout.EdgeColor(:))) / (numberOfLimitsCalculationsSoFar+1);
                
            end
            
            numberOfLimitsCalculationsSoFar = numberOfLimitsCalculationsSoFar +1;
        end

        numberOfRunsSofar = numberOfRunsSofar + 1;
    end
    
    % plot spectra
    if PLOTSPECTRA
        for ch=1:length(PowerChannelsIdx)
            pwr(ch,:) = squeeze(EEG.CAT.Conn.S(PowerChannelsIdx(ch),PowerChannelsIdx(ch),:,:));
        end
        offset = 10*cumsum(ones(length(PowerChannelsIdx),1));
        h=plot(haxpwr,bsxfun(@plus,pwr,offset)');
        set(h,'xdata',EEG.CAT.Conn.freqs);
        axis(haxpwr,'tight');
%         drawnow;
    end
    
    if BENCHMARK
        fprintf('model fitting: %0.5g\n',tfit);
        fprintf('rendering    : %0.5g\n',tdraw);
        fprintf('total        : %0.5g\n',tfit+tdraw);
    end
    
    
    if CHECK_WHITENESS
        [whitestats PC stability] = pop_est_validateMVAR(EEG,0,'checkWhiteness',true,'whitenessCriteria',{'ACF'},...
            'checkConsistency',false,'checkStability',true,'alpha',0.05,'prctWinToSample',100,'verb',0,'plot',0);
        fprintf('\n\n------------------------------\n');
        fprintf('whiteness test:\n');
        fprintf('ACF: p=%0.10g (%s)\n',whitestats{1}.acf.pval,fastif(whitestats{1}.acf.w,'white','not white'));
        %         fprintf('LBP: p=%0.10g (%s)\n',whitestats{1}.ljungbox.pval,fastif(whitestats{1}.ljungbox.w,'white','not white'));
        fprintf('Stability: lambda=%0.10g (%s)\n',max(stability{1}.lambda),fastif(stability{1}.stability,'stable','not stable'));
        fprintf('------------------------------\n\n');
    end
    
    figh = handles.figurehandle;
    
    bmvars = tmp.BMopts.vars;
    mode = 'render';
    
    
    if numberOfRunsSofar == 1
        set(figh,'CloseRequestFcn','assignin(''base'',''closed'',true); delete(gcbf);', ...
            'MenuBar','figure', 'ToolBar','none', 'Name','BrainMovie3D');
    end
    
end

delete(timerfind('Tag','StreamViewerTimer'));
close all force














%% TODO:

% ___ compute channel correlations (offline) and use to interpolate bad
% channels (flt_repair_channels). Can use the dataset in
% /Online/MikeChi/lastdata.set to pre-compute correlations
% ___ band-pass filter data
% ___


