%% RPDC example from (Ex 3.2, eq 16-20) Schelter B, Timmer J, Eichler M (2009) Assessing the strength of directed influences among neural signals using renormalized partial directed coherence. Journal of neuroscience methods 179:121-30 Available at: http://www.ncbi.nlm.nih.gov/pubmed/19428518.

SamplingRate = 100; % Hz

Nl = 3000;               % length of the dataset (5 minutes)
Nr = 1;                  % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
TrueModelOrder = 2;      % VAR model order

expr = { ...
    'x1(t) = 1.9*x1(t-1) + -0.999*x1(t-2) + e1(t)' ...
    'x2(t) = 0.9*x2(t-2) + -0.2*x1(t-1)   + e2(t)' ...
    'x3(t) = -0.3*x3(t-1)+ 0.4*x4(t-1)    + -0.3*x5(t-2) + e3(t)' ...
    'x4(t) = 1.3*x4(t-1) + -0.7*x4(t-2)   + e4(t)' ...
    'x5(t) = 0.7*x5(t-2) + 0.3*x1(t-1)    + e5(t)' ...
    };

% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,TrueModelOrder);


%% STEP 2: Simulate the VAR process

[A] = sim_genTVARcoeffs(Aproto,TrueModelOrder,Nl,'NumSamplesToDiscard',ndisc,'Verbose',true);

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
EEG.setname = 'Sim';
EEG.condition = 'Sim';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pre-process
EEG = pop_pre_prepData(EEG,'nogui','SignalType',{'Channels' 'ConvertChanlocs2Dipfit' []}, 'NormalizeData',[]);

%% Fit the model
WindowLengthSec = EEG.xmax;
WindowStepSec   = 1;
ModelOrder = 10;

EEG = pop_est_fitMVAR(EEG,'nogui','algorithm', 'ARfit','morder',ModelOrder,'winlen',WindowLengthSec,'winstep',WindowStepSec,'epochTimeLims',[],'prctWinToSample',100,'normalize',[],'detrend',[],'timer',0,'verb',2);

%% Create "Ground Truth" model
EEGtrue = EEG;
EEGtrue.CAT.MODEL.AR = A(1);
EEGtrue.CAT.MODEL.PE = {C};
EEGtrue.condition    = 'Truth';
EEGtrue.setname      = 'Truth';

%% Estimate connectivity
ConnMethods = {'S','nPDC','RPDC','PDC'};
Freqs       = linspace(0,EEG.srate/2,50);

EEG     = pop_est_mvarConnectivity(EEG,'nogui','connmethods',ConnMethods,'verb',true,'freqs',Freqs,'ConvertSpectrumToDecibels',true);
EEGtrue = pop_est_mvarConnectivity(EEGtrue,'nogui','connmethods',ConnMethods,'verb',true,'freqs',Freqs,'ConvertSpectrumToDecibels',true);
