% load dataset
% rawdata = io_loadset('eb72_continuous.set'); %io_loadset('bcilab:/userdata/tutorial/flanker_task/12-08-001_ERN.vhdr');

% template EEG set with ICA weights, etc
% EEG = pop_loadset('eb72_continuous.set');   % '/Users/timmullen/Documents/WORK/SIFT/Data/eb79RespWrong.set'

% [1 3 4 6 8 11 16]

% % read from the MINDO device...
run_readmindo('mystream');

%%
Components = [1 3 4 8 16]; %[2:7]; %[1 6 8 11 16]; %[8 11 13 19 20 23 38 39];   % these are the components/channels to which we'll fit our multivariate model
WindowLengthSec = 1;                   % sliding window length in seconds
WindowStepSizeSec = WindowLengthSec;                % sliding window step size in seconds
NewSamplingRate = [];                    % new sampling rate (if downsampling)
% EpochTimeRange = [-1 1.25];              % this is the time range (in seconds) to analyze (relative to event at t=0)
ComponentsToKeep = strtrim(cellstr(num2str(Components')));  % convert list of components to cell array of strings
ModelOrder = 10;

pwrch = 5;

% s.dipfit      = EEG.dipfit;
% s.chanlocs    = EEG.chanlocs;
% s.icaweights  = EEG.icaweights*EEG.icasphere;
% s.icasphere   = eye(length(Components));
% s.icawinv     = EEG.icawinv;
% s.icachansind = EEG.icachansind;
% s.chaninfo    = EEG.chaninfo;
% 
% % delete components
% s.icaweights(setdiff_bc(1:EEG.nbchan,Components),:) = [];
% s.icawinv(:,setdiff_bc(1:EEG.nbchan,Components)) = [];
% 
% clear EEG;

closed = false;
hdata=figure('DeleteFcn',@(varargin)evalin('base','closed=true;'));
haxsrc = axes;

hsrcdata =figure('DeleteFcn',@(varargin)evalin('base','closed=true;'));
haxdata = axes;

hpower =figure('DeleteFcn',@(varargin)evalin('base','closed=true;'));
haxpwr = axes;


%%
% play it back in the background (updated at 4 Hz)
% run_readdataset('mystream',rawdata,25);
% or e.g.: run_readbiosemi('UpdateFrequency',4,'SamplingRate',256);





connmethods = {'dDTF08','S'};

BMfg = struct([]);

BENCHMARK = true;
CHECK_WHITENESS = true;


init = true;
mode = 'init_and_render';

% let the buffer accumulate a bit
pause(1)

EEG = onl_peek('mystream',1);


% poll the buffer periodically
while 1
    
    
    % grab the last 500 ms of data from the stream
    EEG = onl_peek('mystream',WindowLengthSec);
    

    EEG = onl_init_SIFT(EEG,Components,ComponentsToKeep,true);
    
    tic
    
    % fit MVAR model
    EEG = pop_est_fitMVAR(EEG,0,'normalize',{'time'},'algorithm','gladmm','morder',ModelOrder,'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'epochTimeLims',[EEG.xmin EEG.xmax],'verb',0,'gladmm',struct('lambda',0.0005,'verb',0));
    
    % get connectivity
    EEG.CAT.Conn = est_mvarConnectivity(EEG,EEG.CAT.MODEL,'connmethods',connmethods,'freqs',1:40,'verb',0);
    
    if isfield(EEG.CAT.Conn,'S')
        EEG.CAT.Conn = hlp_absvalsq(EEG.CAT.Conn,{'S'});
        EEG.CAT.Conn.S = 10*log10(EEG.CAT.Conn.S);
    end
    
    tfit = toc;
    
    if init
        % initialize brainmovie
        hbmcfg = gui_causalBrainMovie3D_online(EEG,EEG.CAT.Conn,struct('arg_direct',0),BMCFG);
        waitfor(hbmcfg,'UserData','init');
        %         BMCFG.BMopts.speedy = true;
        bmoptschanged = true;
        init = false;
        figh = [];
        bmvars = [];
        
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
    
    
    tic
    
    opts = hlp_mergeVarargin(BMCFG.BMopts,'figurehandle',figh,'mode',mode,'speedy',true,'vars',bmvars);
    [tmp handles] = vis_causalBrainMovie3D(EEG,EEG.CAT.Conn,BMCFG,'timeRange',[],'BMopts',opts);
    
    tdraw = toc;
    
    % plot srcdata
    plot(haxsrc,bsxfun(@plus,EEG.CAT.srcdata',3000*cumsum(ones(length(Components),1))')); drawnow;
    
    % plot original channel data
    plot(haxdata,bsxfun(@plus,EEG.data',3000*cumsum(ones(EEG.nbchan,1))')); drawnow;

    % plot spectrum
    plot(haxpwr,EEG.CAT.Conn.freqs,squeeze(EEG.CAT.Conn.S(pwrch,pwrch,:,:))); drawnow;
    
    
    figh = handles.figurehandle;
    bmvars = tmp.BMopts.vars;
    mode = 'render';
    
    
    if BENCHMARK
        fprintf('model fitting: %0.5g\n',tfit);
        fprintf('rendering    : %0.5g\n',tdraw);
        fprintf('total        : %0.5g\n',tfit+tdraw);
    end
    
        
    if CHECK_WHITENESS
        [whitestats PC stability] = pop_est_validateMVAR(EEG,0,'checkWhiteness',true,'whitenessCriteria',{'Ljung-Box','ACF'},...
                                                'checkConsistency',false,'checkStability',false,'alpha',0.05,'prctWinToSample',100,'verb',0,'plot',0);
        fprintf('\n\n------------------------------\n');
        fprintf('whiteness test:\n');
        fprintf('ACF: p=%0.10g (%s)\n',whitestats{1}.acf.pval,fastif(whitestats{1}.acf.w,'white','not white'));
        fprintf('LBP: p=%0.10g (%s)\n',whitestats{1}.ljungbox.pval,fastif(whitestats{1}.ljungbox.w,'white','not white'));
        fprintf('------------------------------\n\n');
    end
    
    
%     pause(0.01);
end


% onl_clear

