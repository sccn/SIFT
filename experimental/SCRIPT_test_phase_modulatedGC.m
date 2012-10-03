%% PHASE-MODULATED GRANGER CAUSALITY SIMULATION

%% Example 1b: trivariate coupled oscillators with non-stationary (1-Hz sinusoidal) coupling dynamics

%            X1
% (CONST)   /  \    (7 HZ)
%          V    V
%         X2    X3

SamplingRate = 100;      % Hz
Nl = 500;                % length of each epoch (samples)
Nr = 50;                % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
ModelOrder = 2;          % VAR model order
% ModelOrder = 4;        % Model Order we will use


% this is for the est_aamp approach
% TODO need to make x2->x3 be modulated by a theta-oscillatory AR process
% x4 which is also in the system

f0 = 30;                 % central oscillation frequency
expr = {...
    ['x1(t) = ' sim_dampedOscillator(f0,10,100,1) '                                    + e1(t)'] ... 
    ['x2(t) = ' sim_dampedOscillator(f0,2,100,2) ' + -0.1*x1(t-2)                      + e2(t)'] ...
    ['x3(t) = ' sim_dampedOscillator(f0,2,100,3) ' + {fastif(t<250+1000,0,1)*(0.3*sin(2*pi*7*t/100)+0.3+0*randn(1))}*x1(t-2) + e3(t)'] ...
};

% this is for the 
% f0 = 10;              % central oscillation frequency
% expr = {...
%     ['x1(t) = ' sim_dampedOscillator(f0,10,100,1) '                                    + e1(t)'] ... 
%     ['x2(t) = ' sim_dampedOscillator(f0,2,100,2) ' + -0.1*x1(t-2)                      + e2(t)'] ...
%     ['x3(t) = ' sim_dampedOscillator(f0,2,100,3) ' + {0.3*sin(2*pi*fastif(t<250+1000,7,4)*t/100)+0.3+0.1*randn(1)}*x1(t-2) + e3(t)'] ...
% };

% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,ModelOrder);



%% STEP 2: Simulate the VAR process

[A] = sim_genTVARcoeffs(Aproto,ModelOrder,Nl,'NumSamplesToDiscard',ndisc,'Verbose',true);

%% STEP 3: generate data from the VAR model

% Specify the noise covariance matrix. 
% Sigma is the noise variance.
sigma = 1;
M = size(A{1},1);
C = sigma*eye(M);             

% % hyperbolic secant noise
% data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,2/pi,0,'hsec'),[2 1 3]);

% laplacian noise (generalized gaussian)
% data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,1,1,'gengauss'),[2 1 3]);

% gaussian noise
data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,2,1,'gengauss'),[2 1 3]);


%% STEP 4: Create EEGLAB dataset and source potentials

EEG = eeg_emptyset;
EEG.data = data;
[EEG.icaweights EEG.icasphere EEG.icawinv] = deal(eye(M));
EEG.icaact = [];
EEG.srate = SamplingRate;
EEG.times = ((0:(Nl-1))/SamplingRate)*1000;   % ms
EEG.pnts = Nl;
EEG.trials = Nr;
EEG.xmin = EEG.times(1);
EEG.xmax = EEG.times(end)/1000;  % sec
EEG.nbchan = M;
EEG.setname = 'VAR Simulation';
EEG.condition = 'VAR Simulation';

EEG = eeg_checkset(EEG);
pop_eegplot(EEG);

ALLEEG = EEG;
CURRENTSET = length(ALLEEG);
% eeglab redraw;

%%

WindowLengthSec = 1;
WindowStepSizeSec = 0.01;
ModelOrder = 2;

[EEG prepcfg] = pre_prepData('ALLEEG',EEG,'VerbosityLevel',2,'NormalizeData',{'Method',{'time'}});

EEG = est_aamp(EEG,'ampband',[f0-0.25*f0 f0+0.25*f0],'normalize',{'method','time'});

[EEG modfitcfg] = pop_est_fitMVAR(EEG,0,'algorithm','vieira-morf','morder',ModelOrder,'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'verb',1);

[EEG conncfg] = pop_est_mvarConnectivity(EEG,'verb',true,'freqs',(1 : 0.5 : 25),'connmethods',{'dDTF08','nPDC','pCoh','S'},'absvalsq',true,'spectraldecibels',true);

%% Visualize Results

[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,'MatrixLayout',{'Partial','UpperTriangle', 'dDTF08', 'LowerTriangle','dDTF08','Diagonal','S'},'SourceMarginPlot','none','ColorLimits',100,'FrequencyScale','linear','Baseline',[],'Smooth2D',false);

%% integrate over frequencies of interest

connmethod = 'dDTF08';

FrequencyRangeInHz = [f0-0.25*f0 f0+0.25*f0];  % integrate around central frequency using a 25% fractional bandwidth

CollapsedConn = hlp_filterConns(EEG.CAT.Conn,'connmethods',{connmethod},'method',{'freq','net'},'frange',FrequencyRangeInHz,'freqdim',3,'timedim',4);
% if ~isempty(TimeRangeInSeconds), EEG2.CAT.Conn.erWinCenterTimes = 1; end`
%     
% if ~isempty(FrequencyRangeInHz), EEG2.CAT.Conn.freqs = 1; end

[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',CollapsedConn, 'SourceMarginPlot','none','MatrixLayout',{'Partial','UpperTriangle', connmethod, 'LowerTriangle',connmethod,'Diagonal','S'},'ColorLimits',100,'FrequencyScale','linear','Smooth2D',false);

%% spectral analysis of granger-causal coupling

% WARNING, time-varying GC is irregularly sampled, thus we need to use FFT
% (or multitaper) with irregular sampling (or correct the sampling via
% interpolation)

[M M F T] = size(CollapsedConn.(connmethod));  % Note: F is singleton
ConnData = reshape(squeeze(CollapsedConn.(connmethod)),M*M,T);  % [chan-pairs x time]
ConnData = ConnData.'; % [time x chan-pairs]

% apply upsampling interpolation (if window step size is irregular)
% ConnData = interp1(CollapsedConn.winCenterTimes,ConnData,linspace(CollapsedConn.winCenterTimes(1),CollapsedConn.winCenterTimes(end),2*length(CollapsedConn.winCenterTimes)));

BW      = 2;  % Hz
WinLen  = 0.5;  % Sec
params = struct( ...
            'tapers',[BW WinLen 1], ...
            'pad', 2, ...
            'Fs', 1/EEG.CAT.MODEL.winstep, ...
            'err',[1 0.05], ...
            'detrend','linear', ...
            'fpass', [0 50] ...
            );
[S,t,f,Serr] = mtspecgramc(ConnData,[WinLen 0.01],params);
S = permute(S,[3 2 1]); % [chan-pairs x freq x time]
PMGC_Conn.PMGC = reshape(S,M,M,length(f),length(t));  % [chans x chans x freq x time]
PMGC_Conn.winCenterTimes = CollapsedConn.winCenterTimes(1)+t;
PMGC_Conn.erWinCenterTimes = CollapsedConn.erWinCenterTimes(1)+t;
PMGC_Conn.freqs = f;
[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',PMGC_Conn, 'SourceMarginPlot','none','MatrixLayout',{'Partial','UpperTriangle', 'PMGC', 'LowerTriangle','PMGC','Diagonal','none'},'ColorLimits',99.97,'FrequencyScale','linear','Smooth2D',false);


%% example for TA447 data
[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',PMGC_Conn, 'SourceMarginPlot','none','MatrixLayout',{'Partial','UpperTriangle', 'PMGC', 'LowerTriangle','PMGC','Diagonal','none'},'ColorLimits',99.7,'FrequencyScale','linear','Smooth2D',false,'TimesToPlot',[-1.5 5],'Baseline', [],'EventMarkers',{{0,'r',':',2} {tM/1000, 'g', '-', 2}}, 'NodeLabels',anatnames);


%% original TA447 data
[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',CollapsedConn,'Stats',[],'TimesToPlot',[-1.5 3],'Baseline', [],'EventMarkers',{{0,'g',':',2} {tM/1000, 'b', '-', 2}}, 'NodeLabels',anatnames,'SourceMarginPlot','none','MatrixLayout',{'Partial','UpperTriangle', 'dDTF08', 'LowerTriangle','dDTF08','Diagonal','none'},'Thresholding',{'Statistics','ThresholdingMethod','pval','PlotConfidenceIntervals',true,'AlphaSignificance',0.05},'ColorLimits',[],'FrequencyScale','linear','Smooth2D',false);


% newtimef(data, frames, tlimits, SamplingRate, cycles,'key1',value1, 'key2',value2);


