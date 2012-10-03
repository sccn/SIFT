%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to simulate time-varying MVAR process                        %%%
%%% Author: Tim Mullen, 2011, SCCN, INC, UCSD                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% STEP 1: Define the parameters for the system of coupled oscillators
%  This simulation creates a simualted "seizure" with time-varying coupling
%  between clusters of sources which switches between 4 different stages.

%% Simple coupled oscillator example

SamplingRate = 1; % Hz
Nl = 3000;
Nr = 1;                  % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
TrueModelOrder = 2;          % VAR model order
ModelOrder = 2; 50; %200; 
 
expr = {...
    'x1(t) = 0.6*x1(t-1) +  0.65*x2(t-2)+ e1(t)' ... 
    'x2(t) = 0.5*x2(t-1) + -0.3*x2(t-2) + e2(t)' ...
};

% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,TrueModelOrder);


%% RPDC example from (Ex 3.1, eq 11-15) Schelter B, Timmer J, Eichler M (2009) Assessing the strength of directed influences among neural signals using renormalized partial directed coherence. Journal of neuroscience methods 179:121-30 Available at: http://www.ncbi.nlm.nih.gov/pubmed/19428518.

% SamplingRate = 1; % Hz
% 
% Nl = 3000;  % length of the dataset (5 minutes)
% Nr = 1;                  % number of trials (realisations)
% ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
% TrueModelOrder = 3;          % VAR model order
% ModelOrder = 50;
%
% % write the system of equations for coupled damped oscillators
% expr = { ...
%     'x1(t) = 0.9*x1(t-1)  + 0.3*x2(t-2)  + e1(t)' ...
%     'x2(t) = 1.3*x2(t-1)  + -0.8*x2(t-2) + e2(t)' ...
%     'x3(t) = 0.3*x1(t-2)  + 0.6*x2(t-1)  + e3(t)' ...
%     'x4(t) = -0.7*x4(t-3) + -0.7*x1(t-3) + 0.3*x5(t-3) + e4(t)' ...
%     'x5(t) = 1*x5(t-1)    + -0.4*x5(t-2) + 0.3*x4(t-2) + e5(t)' ...
%     };
%   
% % create prototype VAR structure
% Aproto = sim_genVARModelFromEq(expr,TrueModelOrder);


%% RPDC example from (Ex 3.2, eq 16-20) Schelter B, Timmer J, Eichler M (2009) Assessing the strength of directed influences among neural signals using renormalized partial directed coherence. Journal of neuroscience methods 179:121-30 Available at: http://www.ncbi.nlm.nih.gov/pubmed/19428518.

SamplingRate = 1; % Hz

Nl = 3000;  % length of the dataset (5 minutes)
Nr = 500;                  % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
TrueModelOrder = 2;          % VAR model order
ModelOrder = 10; 

expr = { ...
    'x1(t) = 1.9*x1(t-1) + -0.999*x1(t-2) + e1(t)' ...
    'x2(t) = 0.9*x2(t-2) + -0.2*x1(t-1)   + e2(t)' ...
    'x3(t) = -0.3*x3(t-1)+ 0.4*x4(t-1)    + -0.3*x5(t-2) + e3(t)' ...
    'x4(t) = 1.3*x4(t-1) + -0.7*x4(t-2)   + e4(t)' ...
    'x5(t) = 0.7*x5(t-2) + 0.3*x1(t-1)    + e5(t)' ...
    };

% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,TrueModelOrder);


%% PDC Example (eq 5) from Schelter B, Winterhalder M, Eichler M, Peifer M, Hellwig B, Guschlbauer B, Lucking CH, Dahlhaus R, Timmer J (2005) Testing for directed influences among neural signals using partial directed coherence. Journal of neuroscience methods 152:210-9 Available at: http://www.ncbi.nlm.nih.gov/pubmed/16269188.

SamplingRate = 1; % Hz
Nl = 500;
Nr = 100;                  % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
TrueModelOrder = 4;          % VAR model order
ModelOrder = 4; %200; 
 
expr = {...
    'x1(t) = 0.6*x1(t-1) +  0.65*x2(t-2)+e1(t)' ... 
    'x2(t) = 0.5*x2(t-1) + -0.3*x2(t-2)+ -0.3*x3(t-4) +0.6*x4(t-1) + e2(t)' ...
    'x3(t) = 0.8*x3(t-1) + -0.7*x3(t-2)+ -0.1*x5(t-3) + e3(t)' ...
    'x4(t) = 0.5*x4(t-1) +  0.9*x3(t-2)+ 0.4*x5(t-2)  + e4(t)' ...
    'x5(t) = 0.7*x5(t-1) + -0.5*x5(t-2)+ -0.2*x3(t-1) + e5(t)' ...
};

% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,TrueModelOrder);


%% STEP 2: 

[A] = sim_genTVARcoeffs(Aproto,ModelOrder,Nl,'NumSamplesToDiscard',ndisc,'Verbose',true);

%% STEP 3: generate data from the VAR model

M = size(A{1},1);
C = 1*eye(M);             % specify the noise covariance matrix

% % hyperbolic secant noise
% data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,2/pi,0,'hsec'),[2 1 3]);

% laplacian noise (generalized gaussian)
% data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,1,1,'gengauss'),[2 1 3]);

% gaussian noise
data = permute(tvarsim(10*ones(1,M),A,C,[Nl Nr],ndisc,2,1,'gengauss'),[2 1 3]);

%%
data =  permute(arsim(100*ones(1,M),A{1},C,[Nl 1],0),[2 1 3]);

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

EEG = eeg_checkset(EEG);
% pop_eegplot(EEG);


%% 
[EEG prepcfg] = pre_prepData('ALLEEG',EEG,'VerbosityLevel',2,'NormalizeData','off');  %,'NormalizeData','off'

%%
WindowLengthSec = Nl/SamplingRate;
WindowStepSizeSec = Nl/SamplingRate;

[EEG modfitcfg] = pop_est_fitMVAR(EEG,0,'algorithm','vieira-morf','morder',ModelOrder,'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'verb',1);

%%

freqs = linspace(0.01,0.5-0.01,100);
[EEG conncfg] = pop_est_mvarConnectivity(EEG,'verb',false,'freqs',freqs,'connmethods',{'nPDC','RPDC'},'absvalsq',true,'spectraldecibels',false);  %'RPDC','S','nPDC','nDTF'

%% Stats

% surrogate distributions
PConn_phase=stat_surrogate('ALLEEG',EEG,'configs',struct('prepData',prepcfg,'fitMVAR',modfitcfg,'mvarConnectivity',conncfg),'mode',{'PhaseRand', 'nperms', 500},'AutoSave',{'savefname','SIFT_bootstrap','savefrq',10},'connmethods',{'RPDC'},'verb',2);

PConn_boot = stat_surrogate('ALLEEG',EEG,'configs',struct('prepData',prepcfg,'fitMVAR',modfitcfg,'mvarConnectivity',conncfg),'mode',{'Bootstrap', 'nperms', 500},'AutoSave',{'savefname','SIFT_bootstrap','savefrq',10},'connmethods',{'RPDC','nPDC'},'verb',2);

PConn_invjack = stat_surrogate('ALLEEG',EEG,'configs',struct('prepData',prepcfg,'fitMVAR',modfitcfg,'mvarConnectivity',conncfg),'mode',{'InverseJacknife'},'AutoSave',{'savefname','SIFT_bootstrap','savefrq',10},'connmethods',{'RPDC','nPDC','nDTF'},'verb',2);

PConn_jack = stat_surrogate('ALLEEG',EEG,'configs',struct('prepData',prepcfg,'fitMVAR',modfitcfg,'mvarConnectivity',conncfg),'mode',{'Jacknife'},'AutoSave',{'savefname','SIFT_bootstrap','savefrq',10},'connmethods',{'RPDC','nPDC','nDTF'},'verb',2);


% compute stats (REPLACE WITH stat_surrogateStats)

Stats = stat_bootSigTest('PConn',EEG.CAT.Conn,'null',PConn_phase,'ConnectivityMethods',{'RPDC'},'mcorrection','fdr','statcondargs',{'tail','one'});

Stats = stat_bootSigTest('PConn',EEG.CAT.Conn,'null',PConn,'ConnectivityMethods',{'RPDC'},'mcorrection','fdr','statcondargs',{'tail','one'});


%%
% compute analytic stats
Stats = stat_analyticStats('ALLEEG',EEG,'Estimator',{'RPDC','nPDC'},'Alpha', 0.05,'verb',false);

%%

[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,'Stats',Stats,'MatrixLayout',{'Partial','UpperTriangle', 'RPDC', 'LowerTriangle','RPDC','Diagonal','none'},'plotci',true,'ColorLimits',[0 0.025],'Baseline',[],'Smooth2D',false,'Thresholding',{'Statistics','AlphaSignificance',0.05},'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[1 1 1]);

%% no stats

vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,'MatrixLayout',{'Partial','UpperTriangle', 'RPDC', 'LowerTriangle','RPDC','Diagonal','none'},'plotci',false,'ColorLimits',100,'Baseline',[],'Smooth2D',false,'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[1 1 1]);


%% no stats (multiple conditions)

vis_TimeFreqGrid('ALLEEG',EEG,'Conn',[EEG(2).CAT.Conn EEG(1).CAT.Conn],'MatrixLayout',{'Partial','UpperTriangle', 'nPDC', 'LowerTriangle','nPDC','Diagonal','none'},'plotci',false,'ColorLimits',100,'Baseline',[],'Smooth2D',false,'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[1 1 1]);

%%

vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,'Stats',Stats,'MatrixLayout',{'Partial','UpperTriangle', 'nPDC', 'LowerTriangle','nPDC','Diagonal','none'},'plotci',false,'ColorLimits',100,'Baseline',[],'Smooth2D',false,'Thresholding',{'Statistics','ThresholdingMethod','pval','AlphaSignificance',0.05},'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[1 1 1]);

%% Test new plotting routines

vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,'Stats',Stats,'MatrixLayout',{'Partial','UpperTriangle', 'nPDC', 'LowerTriangle','nPDC','Diagonal','none'},'plotci',false,'ColorLimits',100,'Baseline',[],'Smooth2D',false,'Thresholding',{'Statistics','ThresholdingMethod','pval','AlphaSignificance',0.05},'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[1 1 1],'PatchColor','y');
vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,'Stats',Stats,'MatrixLayout',{'Partial','UpperTriangle', 'nPDC', 'LowerTriangle','nPDC','Diagonal','none'},'plotci',false,'ColorLimits',100,'Baseline',[],'Smooth2D',false,'Thresholding',{'Statistics','ThresholdingMethod','thresh','AlphaSignificance',0.05},'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[1 1 1],'PatchColor','y');

%%

[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,'Stats',Stats,'MatrixLayout',{'Partial','UpperTriangle', 'nPDC', 'LowerTriangle','nPDC','Diagonal','none'},'plotci',true,'ColorLimits',99.7,'Baseline',[],'Smooth2D',false,'Thresholding',{'Statistics','ThresholdingMethod','thresh','AlphaSignificance',0.05},'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[1 1 1]);

%%

Stats = stat_bootSigTest('PConn',EEG.CAT.Conn,'null',stat_getBaselineDistrib(PConn_boot,[-1 -0.25]),'ConnectivityMethods',{'RPDC'},'mcorrection','fdr','statcondargs',{'tail','one'});