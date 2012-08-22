% load dataset
rawdata = io_loadset('eb72_continuous.edf'); % 'eb72_continuous.edf'  %io_loadset('bcilab:/userdata/tutorial/flanker_task/12-08-001_ERN.vhdr');

% template EEG set with ICA weights, etc
EEG = pop_loadset('eb72_continuous.set');   

% [8 11 12 13 14 15 18 19 20 23 24 28 38 39 60 65]
% [8 11 13 19 20 23 38 39];
Components = [8 11 13 19 20 23 38 39];  %[8 11 13 19 20];   % these are the components/channels to which we'll fit our multivariate model
WindowLengthSec = 0.5;                   % sliding window length in seconds
WindowStepSizeSec = 0.03; WindowLengthSec;                % sliding window step size in seconds
NewSamplingRate = [];                    % new sampling rate (if downsampling)
% EpochTimeRange = [-1 1.25];              % this is the time range (in seconds) to analyze (relative to event at t=0)
ComponentsToKeep = strtrim(cellstr(num2str(Components')));  % convert list of components to cell array of strings
ModelOrder = 10;


s.dipfit      = EEG.dipfit;
s.chanlocs    = EEG.chanlocs;
s.icaweights  = EEG.icaweights*EEG.icasphere;
s.icasphere   = eye(length(Components));
s.icawinv     = EEG.icawinv;
s.icachansind = EEG.icachansind;
s.chaninfo    = EEG.chaninfo;

% delete components
% s.icaweights(setdiff(1:EEG.nbchan,Components),:) = [];
% s.icawinv(:,setdiff(1:EEG.nbchan,Components)) = [];

clear EEG;



%%
% play it back in the background (updated at 4 Hz)
run_readdataset('mystream',rawdata,25);


% read from the MINDO device...
% run_readmindo('mystream');



% or e.g.: run_readbiosemi('UpdateFrequency',4,'SamplingRate',256);

% brainmovie figure
% figh = figure;

% BMcfg = struct(...
%         'connmethod','nPDC' , ...
%         'freqsToCollapse',[3:10] , ...
%         'collapsefun','integrate', ...
%         'resample',0, ...
%         'subtractconds',0, ...
%         'nodelabels',ComponentsToKeep, ...
%         'nodesToExclude',[], ...
%         'edgeColorMapping','Connectivity', ...
%         'edgeSizeMapping','ConnMagnitude', ...
%         'nodeColorMapping','Outflow', ...
%         'nodeSizeMapping','Outflow', ...
%         'baseline',[], ...
%         'normalize',1, ...
%         'prcthresh',[], ...
%         'absthresh',[], ...
%         'footerPanelSpec','off', ...
%           'BMopts', struct( ...
%               'visible', 'on',    ...
%               'latency', [],      ...
%               'frames', [],       ...
%               'figurehandle', figh, ...
%               'speedy', true,     ...
%               'rotationpath3d', 'none', ...
%               'project3d', 'off', ...
%               'plotCortex', struct( ...
%                             'cortexVolumeFile', 'standard_BEM_vol.mat', ...
%                             'cortexTransparency', 0.7 ...
%                             ), ...
%               'opengl', 'on', ...
%               'flashes', [], ...
%               'square', 'on', ...
%               'caption', 'off', ...
%               'showLatency', 1, ...
%               'dispRT', 0, ...
%               'backcolor', [0 0 0], ...
%               'graphColorAndScaling', struct(...
%                         'nodeSizeLimits', [0.1 1],  ...
%                         'nodeColorLimits', [0 1],   ...
%                         'edgeSizeLimits', [0.1 1],  ...
%                         'edgeColorLimits', [0 1],   ...
%                         'nodeSizeDataRange', [],    ...
%                         'nodeColorDataRange', [],   ...
%                         'edgeSizeDataRange', [],    ...
%                         'edgeColorDataRange', [],   ...
%                         'centerDataRange', 0,       ...
%                         'edgeColormap', jet(64), ...
%                         'nodeColormap', jet(64), ...
%                         'diskscale', 0.3, ...
%                         'magnify', 1 ...
%                         ), ...
%                 'outputFormat', struct(...
%                                 'framefolder', '', ...
%                                 'framesout', 'jpg', ...
%                                 'moviename', '', ...
%                                 'movieopts', {'videoname' ''}, ...
%                                 'size', [600 600] ...
%                                 ), ...
%                 'mri', '', ...
%                 'coordformat', 'spherical', ...
%                 'dipplotopt', {}, ...
%                 'renderBrainMovie', 1, ...
%                 'initonly', true ...
%                 ) ...
%             );

connmethods = {'dDTF08'};

BMfg = struct([]);

BENCHMARK = true;
CHECK_WHITENESS = true;
CROSS_FADE = false;

init = true;
mode = 'init_and_render';

% let the buffer accumulate a bit
pause(1)

EEG = onl_peek('mystream',1);

% datafigure = figure;
% dataax = axes;

% ff = figure;
% ax = axes;


% poll the buffer periodically
while 1
    
    
    % grab the last 500 ms of data from the stream
    EEG = onl_peek('mystream',WindowLengthSec);
    
    
    % insert missing fields (replace with initialization in run_readbiosemi)
    f = fieldnames(s);
    for i = 1:length(f)
        EEG.(f{i}) = s.(f{i});
    end
    
    % pre-process in SIFT
    %     EEG = pre_prepData('ALLEEG',EEG,'VerbosityLevel',0,'NewSamplingRate',NewSamplingRate,'NormalizeData',{'Method',{'time'}},'SelectComponents',{'ComponentsToKeep',ComponentsToKeep});
    %     EEG.times = linspace(EEG.xmin*1000,EEG.xmax*1000,EEG.pnts);
    
    EEG = onl_init_SIFT(EEG,Components,ComponentsToKeep);
    
    tic
    
    % fit MVAR model
    EEG = pop_est_fitMVAR(EEG,0,'normalize',{'time'},'algorithm','gladmm','morder',ModelOrder,'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'epochTimeLims',[EEG.xmin EEG.xmax],'verb',0,'gladmm',struct('lambda',0.0005,'verb',0));
    
    % get connectivity
    EEG.CAT.Conn = est_mvarConnectivity(EEG,EEG.CAT.MODEL,'connmethods',connmethods,'freqs',1:40,'verb',0);
   
    tfit = toc;
    
    
    if CHECK_WHITENESS
        [whitestats PC stability] = pop_est_validateMVAR(EEG,0,'checkWhiteness',true,'whitenessCriteria',{'Ljung-Box','ACF'},...
                                                'checkConsistency',false,'checkStability',false,'alpha',0.05,'prctWinToSample',100,'verb',0,'plot',0);
        fprintf('\n\n------------------------------\n');
        fprintf('whiteness test:\n');
        fprintf('ACF: p=%0.10g (%s)\n',whitestats{1}.acf.pval,fastif(whitestats{1}.acf.w,'white','not white'));
        fprintf('LBP: p=%0.10g (%s)\n',whitestats{1}.ljungbox.pval,fastif(whitestats{1}.ljungbox.w,'white','not white'));
        fprintf('------------------------------\n\n');
    end
    
    if init
        % initialize brainmovie
        hbmcfg = gui_causalBrainMovie3D_online(EEG,EEG.CAT.Conn,struct('arg_direct',0));
        waitfor(hbmcfg,'UserData','init');
        %         BMCFG.BMopts.speedy = true;
        bmoptschanged = true;
        init = false;
        figh = [];
        bmvars = [];
        
        if CROSS_FADE
            % cross-fade feature
            fadetimer = timer;
        end

        
%         eeglines=plot(dataax,EEG.times,squeeze(EEG.CAT.srcdata(1,:,:)));
%         
%         eegplot(EEG.CAT.srcdata,'srate',EEG.srate);
%         eeglines = get(findobj(gcf,'tag','eegaxis'),'children');
        
%         OldConn = EEG.CAT.Conn;
    else
%         figh = BMCFG.BMopts.figurehandle;
        
%         norm(OldConn.dDTF08(:)-EEG.CAT.Conn.dDTF08(:));
%         OldConn = EEG.CAT.Conn;
    end
    
    if bmoptschanged
        mode = 'init_and_render';
        bmoptschanged = false;
        bmvars = [];
%         BMCFG
        
    end
    
    %     fprintf('data norm: %0.5g\n', norm(EEG.data(:)));
    
    %     norm(EEG.CAT.Conn.dDTF08(:))
    
%     plot(ax,EEG.CAT.Conn.freqs, squeeze(EEG.CAT.Conn.dDTF08(1,3,:,1)));
    
    
    % render brainmovie
     
    tic
    
    opts = hlp_mergeVarargin(BMCFG.BMopts,'figurehandle',figh,'mode',mode,'speedy',true,'vars',bmvars);
    [tmp handles] = vis_causalBrainMovie3D(EEG,EEG.CAT.Conn,BMCFG,'timeRange',[],'BMopts',opts);
    
    tdraw = toc;
    
    % testing cross-fade feature
    
    if ~isempty(figh) && CROSS_FADE
%         surfaces = findobj(get(findobj(figh,'tag','brain1'),'children'),'type','surface'); % all surfaces in brain

        surfaces = findobj(figh,'tag','tmpmov','type','surface'); % all edges/nodes
        
        stop(fadetimer);
        set(fadetimer,'Period',round(1000*((tdraw+tfit)/10))/1000, ...  % round(1000*((tdraw+tfit)/10))/1000
            'TimerFcn',{@(obj,evnt,handles,factor) hlp_scaleAlpha(handles,factor),surfaces,0.8}, ...
            'ExecutionMode','fixedDelay','TasksToExecute',15);
        start(fadetimer);
        
%         hlp_scaleAlpha(surfaces,0.5);
    end
        

%     set(eeglines,'ydata',EEG.CAT.srcdata(1,:,:));
    
%     EEG.data = EEG.CAT.srcdata;
%     vis_hist(EEG);

%     eegplot(EEG.CAT.srcdata,'srate',EEG.srate);
    
%     for i=1:length(eeglines)
%         set(eeglines(i),'ydata',EEG.CAT.srcdata(i,:,:));
%     end
    
    figh = handles.figurehandle;
        
    bmvars = tmp.BMopts.vars;
    
    mode = 'render';
    
    
    if BENCHMARK
        fprintf('model fitting: %0.5g\n',tfit);
        fprintf('rendering    : %0.5g\n',tdraw);
        fprintf('total        : %0.5g\n',tfit+tdraw);
    end
%     pause(0.01);
end


% onl_clear

