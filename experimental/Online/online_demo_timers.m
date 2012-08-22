% load dataset
rawdata = io_loadset('eb72_continuous.edf'); % 'eb72_continuous.edf'  %io_loadset('bcilab:/userdata/tutorial/flanker_task/12-08-001_ERN.vhdr');

% template EEG set with ICA weights, etc
EEG = pop_loadset('eb72_continuous.set');   

% [8 11 12 13 14 15 18 19 20 23 24 28 38 39 60 65]
% [8 11 13 19 20 23 38 39];
Components = [8 11 13 19];   % these are the components/channels to which we'll fit our multivariate model
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
run_readdataset('mystream',rawdata,10);
% or e.g.: run_readbiosemi('UpdateFrequency',4,'SamplingRate',256);

connmethods = {'dDTF08'};

BMfg = struct([]);

BENCHMARK = true;
CHECK_WHITENESS = true;

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
maintimer = timer('Period',0.25,'TimerFcn','eval(online_inner_loop)','BusyMode','queue','ExecutionMode','fixedDelay','tag','LoopTimer');


start(maintimer);

%%
stop(timerfind('tag','LoopTimer'));
delete(timerfind('tag','LoopTimer'));
clear maintimer


