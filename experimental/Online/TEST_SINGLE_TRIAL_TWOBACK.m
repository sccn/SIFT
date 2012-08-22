% test continuous data 

EEG = pop_loadset('eb72_continuous.set');

%% 
% [8 11 12 13 14 15 18 19 20 23 24 28 38 39 60 65]
% [8 11 13 19 20 23 38 39]
Components = [8 11]; % these are the components/channels to which we'll fit our multivariate model
WindowLengthSec = 1;                 % sliding window length in seconds
WindowStepSizeSec = 0.1;              % sliding window step size in seconds
NewSamplingRate = [];                  % new sampling rate (if downsampling)
EpochTimeRange = [0 30];              % this is the time range (in seconds) to analyze (relative to event at t=0)
ComponentsToKeep = strtrim(cellstr(num2str(Components')));  % convert list of components to cell array of strings
ModelOrder = 12;

[EEG prepcfg] = pre_prepData('ALLEEG',EEG,'verb',2,'backupOriginalData',0, ...
                             'selectComponents',{'verb' 0 'ComponentsToKeep' ComponentsToKeep}, ...
                             'newtlims',EpochTimeRange, ...
                             'resample',[], ...
                             'filter',[], ...
                             'detrend',[], ...   %   {'arg_selection' 0}   % 'verb' 1 'method' {'linear'} 'piecewise' {'seglength' WindowLengthSec 'stepsize' WindowStepSizeSec 'arg_selection' 1} 'plot' 0 'arg_selection' 1
                             'normalize',[], ...     %  'arg_direct' 0 'verb' 1 'method' {'time'} 'arg_selection' 0
                             'newtrials',[], ...
                             'badsegments',[], ...
                             'equalizetrials',0);

EEG.times = linspace(EpochTimeRange(1),EEG.xmax*1000,EEG.pnts);

%%

connmethods = {'dDTF08','nPDC','S'};

% fit MVAR model
EEG = pop_est_fitMVAR(EEG,0,'normalize',{'time'},'algorithm','gladmm','morder',ModelOrder,'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'epochTimeLims',EpochTimeRange,'verb',1,'gladmm',struct('lambda',10,'verb',1,'max_iter',10000));

% get connectivity
EEG.CAT.Conn = est_mvarConnectivity(EEG,EEG.CAT.MODEL,'connmethods',connmethods,'freqs',1:40,'verb',0);
EEG.CAT.Conn = hlp_absvalsq(EEG.CAT.Conn,{'S','nPDC'});
EEG.CAT.Conn.S = 10*log10(EEG.CAT.Conn.S);



