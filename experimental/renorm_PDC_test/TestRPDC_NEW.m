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
ModelOrder = 200;

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

%% Visualization
cfg = {'plotCondDiff',0,'vismode','TimeXFrequency','MatrixLayout',{'arg_selection' 'Full' 'estimator' 'RPDC'},'clim',100,'timeRange',[14.995 14.995] ,'freqValues',[0      1.02041      2.04082      3.06122      4.08163      5.10204      6.12245      7.14286      8.16327      9.18367      10.2041      11.2245      12.2449      13.2653      14.2857      15.3061      16.3265      17.3469      18.3673      19.3878      20.4082      21.4286       22.449      23.4694      24.4898      25.5102      26.5306       27.551      28.5714      29.5918      30.6122      31.6327      32.6531      33.6735      34.6939      35.7143      36.7347      37.7551      38.7755      39.7959      40.8163      41.8367      42.8571      43.8776       44.898      45.9184      46.9388      47.9592      48.9796           50] ,'windows',[],'pcontour',[],'thresholding',{'arg_selection' 'None'},'baseline',[],'fighandles',[],'smooth',0,'xord',[],'yord',[],'plotorder',[],'topoplot','none','customTopoMatrix',[],'dipplot',{'mri' '' 'coordformat' 'mni' 'dipplotopt' {}},'nodelabels',{'1' '2' '3' '4' '5'},'foilines',[],'foilinecolor',[],'events',{{0 'r' ':' 2}},'freqscale','linear','transform','linear','titleString','','titleFontSize',12,'axesFontSize',10,'textColor',[1 1 1] ,'linecolor',[1 1 1] ,'patchcolor',[1 1 1] ,'colormap',[0 0 0.50667;0 0 0.52;0 0 0.53333;0 0 0.54667;0 0 0.56;0 0 0.57333;0 0 0.58667;0 0 0.6;0 0 0.61333;0 0 0.62667;0 0 0.64;0 0 0.65333;0 0 0.66667;0 0 0.68;0 0 0.69333;0 0 0.70667;0 0 0.72;0 0 0.73333;0 0 0.74667;0 0 0.76;0 0 0.77333;0 0 0.78667;0 0 0.8;0 0 0.81333;0 0 0.82667;0 0 0.84;0 0 0.85333;0 0 0.86667;0 0 0.88;0 0 0.89333;0 0 0.90667;0 0 0.92;0 0 0.93333;0 0 0.94667;0 0 0.96;0 0 0.97333;0 0 0.98667;0 0 1;0 0.013333 1;0 0.026667 1;0 0.04 1;0 0.053333 1;0 0.066667 1;0 0.08 1;0 0.093333 1;0 0.10667 1;0 0.12 1;0 0.13333 1;0 0.14667 1;0 0.16 1;0 0.17333 1;0 0.18667 1;0 0.2 1;0 0.21333 1;0 0.22667 1;0 0.24 1;0 0.25333 1;0 0.26667 1;0 0.28 1;0 0.29333 1;0 0.30667 1;0 0.32 1;0 0.33333 1;0 0.34667 1;0 0.36 1;0 0.37333 1;0 0.38667 1;0 0.4 1;0 0.41333 1;0 0.42667 1;0 0.44 1;0 0.45333 1;0 0.46667 1;0 0.48 1;0 0.49333 1;0 0.50667 1;0 0.52 1;0 0.53333 1;0 0.54667 1;0 0.56 1;0 0.57333 1;0 0.58667 1;0 0.6 1;0 0.61333 1;0 0.62667 1;0 0.64 1;0 0.65333 1;0 0.66667 1;0 0.68 1;0 0.69333 1;0 0.70667 1;0 0.72 1;0 0.73333 1;0 0.74667 1;0 0.76 1;0 0.77333 1;0 0.78667 1;0 0.8 1;0 0.81333 1;0 0.82667 1;0 0.84 1;0 0.85333 1;0 0.86667 1;0 0.88 1;0 0.89333 1;0 0.90667 1;0 0.92 1;0 0.93333 1;0 0.94667 1;0 0.96 1;0 0.97333 1;0 0.98667 1;0 1 1;0.013333 1 0.98667;0.026667 1 0.97333;0.04 1 0.96;0.053333 1 0.94667;0.066667 1 0.93333;0.08 1 0.92;0.093333 1 0.90667;0.10667 1 0.89333;0.12 1 0.88;0.13333 1 0.86667;0.14667 1 0.85333;0.16 1 0.84;0.17333 1 0.82667;0.18667 1 0.81333;0.2 1 0.8;0.21333 1 0.78667;0.22667 1 0.77333;0.24 1 0.76;0.25333 1 0.74667;0.26667 1 0.73333;0.28 1 0.72;0.29333 1 0.70667;0.30667 1 0.69333;0.32 1 0.68;0.33333 1 0.66667;0.34667 1 0.65333;0.36 1 0.64;0.37333 1 0.62667;0.38667 1 0.61333;0.4 1 0.6;0.41333 1 0.58667;0.42667 1 0.57333;0.44 1 0.56;0.45333 1 0.54667;0.46667 1 0.53333;0.48 1 0.52;0.49333 1 0.50667;0.50667 1 0.49333;0.52 1 0.48;0.53333 1 0.46667;0.54667 1 0.45333;0.56 1 0.44;0.57333 1 0.42667;0.58667 1 0.41333;0.6 1 0.4;0.61333 1 0.38667;0.62667 1 0.37333;0.64 1 0.36;0.65333 1 0.34667;0.66667 1 0.33333;0.68 1 0.32;0.69333 1 0.30667;0.70667 1 0.29333;0.72 1 0.28;0.73333 1 0.26667;0.74667 1 0.25333;0.76 1 0.24;0.77333 1 0.22667;0.78667 1 0.21333;0.8 1 0.2;0.81333 1 0.18667;0.82667 1 0.17333;0.84 1 0.16;0.85333 1 0.14667;0.86667 1 0.13333;0.88 1 0.12;0.89333 1 0.10667;0.90667 1 0.093333;0.92 1 0.08;0.93333 1 0.066667;0.94667 1 0.053333;0.96 1 0.04;0.97333 1 0.026667;0.98667 1 0.013333;1 1 0;1 0.98667 0;1 0.97333 0;1 0.96 0;1 0.94667 0;1 0.93333 0;1 0.92 0;1 0.90667 0;1 0.89333 0;1 0.88 0;1 0.86667 0;1 0.85333 0;1 0.84 0;1 0.82667 0;1 0.81333 0;1 0.8 0;1 0.78667 0;1 0.77333 0;1 0.76 0;1 0.74667 0;1 0.73333 0;1 0.72 0;1 0.70667 0;1 0.69333 0;1 0.68 0;1 0.66667 0;1 0.65333 0;1 0.64 0;1 0.62667 0;1 0.61333 0;1 0.6 0;1 0.58667 0;1 0.57333 0;1 0.56 0;1 0.54667 0;1 0.53333 0;1 0.52 0;1 0.50667 0;1 0.49333 0;1 0.48 0;1 0.46667 0;1 0.45333 0;1 0.44 0;1 0.42667 0;1 0.41333 0;1 0.4 0;1 0.38667 0;1 0.37333 0;1 0.36 0;1 0.34667 0;1 0.33333 0;1 0.32 0;1 0.30667 0;1 0.29333 0;1 0.28 0;1 0.26667 0;1 0.25333 0;1 0.24 0;1 0.22667 0;1 0.21333 0;1 0.2 0;1 0.18667 0;1 0.17333 0;1 0.16 0;1 0.14667 0;1 0.13333 0;1 0.12 0;1 0.10667 0;1 0.093333 0;1 0.08 0;1 0.066667 0;1 0.053333 0;1 0.04 0;1 0.026667 0;1 0.013333 0;1 0 0;0.98667 0 0;0.97333 0 0;0.96 0 0;0.94667 0 0;0.93333 0 0;0.92 0 0;0.90667 0 0;0.89333 0 0;0.88 0 0;0.86667 0 0;0.85333 0 0;0.84 0 0;0.82667 0 0;0.81333 0 0;0.8 0 0;0.78667 0 0;0.77333 0 0;0.76 0 0;0.74667 0 0;0.73333 0 0;0.72 0 0;0.70667 0 0;0.69333 0 0;0.68 0 0;0.66667 0 0;0.65333 0 0;0.64 0 0;0.62667 0 0;0.61333 0 0;0.6 0 0;0.58667 0 0;0.57333 0 0;0.56 0 0;0.54667 0 0;0.53333 0 0;0.52 0 0;0.50667 0 0],'backgroundColor',[0 0 0] };
vis_TimeFreqGrid('EEG',EEG,'Conn',EEG.CAT.Conn,cfg{:});
