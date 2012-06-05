
function [EEGprep args] = pre_prepData(varargin)
%
% Preprocess EEG dataset(s) for connectivity analysis. See [1] for
% mathematical details on preprocessing steps.
%
%
% Input:
%
%   ALLEEG:         Array of EEGLAB datasets to preprocess.
%
% Optional:         <'Name',value> pairs
%
%     VerbosityLevel:        Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                            Possible values: 0,1,2
%                            Default value  : 0
%                            Input Data Type: real number (double)
%
%     backupOriginalData:    Keep a nonnormalized copy of the data
%                            Input Data Type: boolean
%
%     SelectComponents:      Select components to analyze
%                            Input Data Type: string
%     -----------------
%
%         VerbosityLevel:    Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                            Possible values: 0,1,2
%                            Default value  : 0
%                            Input Data Type: real number (double)
%
%         ComponentsToKeep:  Select components to analyze
%                            This should be a cell array of strings containing the IDs of components you wish to keep
%                            Input Data Type: boolean
%
%     EpochTimeRange:        [Min Max] Epoch time range (sec)
%                            If blank, use all original epoch length
%                            Input Data Type: real number (double)
%
%     NewSamplingRate:       New sampling rate
%                            Data will be down/upsampled using a zero-phase filter (see 'help resample')
%                            Input Data Type: real number (double)
%
%     FilterData:            Filter specifications
%                            Provide [low hi] (Hz) edges of bandpass filter. [0 f] will provide a low-pass filter at frequency
%                            f. [f 0] will produce a high-pass filter at frequency f.
%                            Input Data Type: real number (double)
%
%     DifferenceData:        Differencing options
%                            A difference filter of order k applied to time series X is defined as [for i=1:k, for t=1:N, X(t) =
%                            X(t)-X(t-1)]
%                            Input Data Type: boolean
%     ---------------
%
%         VerbosityLevel:    Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                            Possible values: 0,1,2
%                            Default value  : 1
%                            Input Data Type: real number (double)
%
%         DifferencingOrder: Differencing order
%                            Number of times to difference data
%                            Input Range  : [0  10]
%                            Default value: 1
%                            Input Data Type: real number (double)
%
%     Detrend:               Detrend or center each epoch
%                            Input Data Type: boolean
%     --------
%
%         VerbosityLevel:    Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                            Possible values: 0,1,2
%                            Default value  : 1
%                            Input Data Type: real number (double)
%
%         DetrendingMethod:  Detrending options
%                            Linear: removes the least-squares fit of a straight line from each trial. Constant: removes the
%                            mean from each trial (centering)
%                            Possible values: 'linear','constant'
%                            Default value  : 'linear'
%                            Input Data Type: boolean
%
%     NormalizeData:         Data normalization
%                            Normalize trials across time, ensemble, or both
%                            Input Data Type: boolean
%     --------------
%
%         VerbosityLevel:    Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                            Possible values: 0,1,2
%                            Default value  : 0
%                            Input Data Type: real number (double)
%
%         Method:            Normalize windows across time, ensemble, or both
%                            Possible values: 'ensemble','time'
%                            Default value  : 'ensemble'
%                            Input Data Type: boolean
%
%     TrialSubsetToUse:      Subset of trial indices to use
%                            Input Data Type: real number (double)
%
%     BadDataSegments:       Intervals of bad data
%                            N x 2 matrix of [lo hi] intervals (seconds) of data within ea. trial to set to nan. Currently only
%                            compatible with vierra_morf MVAR method
%                            Input Data Type: real number (double)
%
%     BalanceTrials:         Equalize the number of trials between two conditions
%                            Input Data Type: boolean
%
% Output:
%
%   EEGprep:        Prepocessed EEG structure(s)
%   g:              Argument specification structure.
%
%
% See Also: pop_pre_prepData()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



% extract some defaults
ALLEEG = arg_extract(varargin,'ALLEEG',1);
MyComponentNames = [];

procChans = false;

if ~isempty(ALLEEG)
    
    % make sure we have ICA activations
    if any(cellfun(@isempty,{ALLEEG.icaweights}))
        res = questdlg2('ICA decomposition not available. Process channels instead?','SIFT Data Preprocessing','Yes','No','Yes');
        if strcmpi(res,'Yes')
            
            procChans = true;
            
            for cond=1:length(ALLEEG)
                % insert 'fake' ICA solution (Identity matrices) and copy EEG data into icaact
                [ALLEEG(cond).icaweights ALLEEG(cond).icasphere ALLEEG(cond).icawinv] = deal(eye(size(ALLEEG(cond).data,1)));
                ALLEEG(cond).icaact = ALLEEG(cond).data;
            end
        else
            EEGprep = ALLEEG;
            return;
        end
    end
    
    for cond=1:length(ALLEEG)
        
        
        if ~isfield(ALLEEG(cond),'CAT') || ~isfield(ALLEEG(cond).CAT,'curComps')
            ALLEEG(cond).CAT.curComps = 1:size(ALLEEG(cond).icaweights,1);
        end
        if  ~isfield(ALLEEG(cond),'CAT') || ~isfield(ALLEEG(cond).CAT,'MODEL')
            ALLEEG(cond).CAT.MODEL = [];
        end
    end
    
    
    if isfield(ALLEEG(1).CAT,'curComponentNames') && ~isempty(ALLEEG(1).CAT.curComponentNames)
        MyComponentNames = ALLEEG(1).CAT.curComponentNames;
    else
        MyComponentNames = ALLEEG(1).CAT.curComps;
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    end
    clear ALLEEG
end

verb = arg_extract(varargin,{'verb','VerbosityLevel'},[],0);

g = arg_define([0 Inf], varargin, ...
    arg_norep('ALLEEG',mandatory), ...
    arg({'verb','VerbosityLevel'},0,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
    arg_nogui({'backupOriginalData'},false,[],'Keep a nonnormalized copy of the data'), ...
    arg_sub({'selectComponents','SelectComponents'},{'MyComponentNames',MyComponentNames,'verb',verb},@pre_selectcomps,'Select components to analyze','cat','Data Selection'), ...
    arg({'newtlims','EpochTimeRange'},[],[],'[Min Max] Epoch time range (sec). If blank, use all original epoch length','cat','Data Selection'), ...
    arg({'resample','NewSamplingRate'},[],[],'New sampling rate. Data will be down/upsampled using a zero-phase filter (see ''help resample'')','cat','Filtering'), ...
    arg({'filter','FilterData'},[],[],'Filter specifications. Provide [low hi] (Hz) edges of bandpass filter. [0 f] will provide a low-pass filter at frequency f. [f 0] will produce a high-pass filter at frequency f.','shape','row','cat','Filtering'), ...
    arg_subtoggle({'diff','DifferenceData'},[],@pre_diffData,'Differencing options. A difference filter of order k applied to time series X is defined as [for i=1:k, for t=1:N, X(t) = X(t)-X(t-1)]','cat','Filtering'), ...
    arg_subtoggle({'detrend','Detrend'},[],@pre_detrend,'Detrend or center each epoch','cat','Filtering'), ...
    arg_subtoggle({'normalize','NormalizeData'},{'method', {'time','ensemble'}, 'verb',verb},@pre_normData,'Data normalization. Normalize trials across time, ensemble, or both','cat','Normalization'), ...
    arg_subtoggle({'aamp','AmplitudeEnvelope'},[],@est_aamp,'Compute amplitude envelope. This will add additional ''pseudochannels'' which represent the amplitude envelopes of original channels over a specified frequency band. Envelopes are computed via a zero-phase filter (eegfilt) followed by a hilbert transform.','cat','Filtering'), ...
    arg({'newtrials','TrialSubsetToUse'},[],[],'Subset of trial indices to use','cat','Data Selection'), ...
    arg_nogui({'badsegments','BadDataSegments'},[],[],'Intervals of bad data. N x 2 matrix of [lo hi] intervals (seconds) of data within ea. trial to set to nan. Currently only compatible with vierra_morf MVAR method','cat','Data Selection'), ...
    arg_nogui({'equalizetrials','BalanceTrials'},false,[],'Equalize the number of trials between two conditions','cat','Data Selection'), ...
    arg_norep({'DO_JACKET','UseJacket'},false,[],'Use Jacket (GPU)','cat','Miscellaneous') ...
    );

% set the EEGLAB memory options to SIFT defaults
if g.verb, fprintf('Setting EEGLAB options to SIFT defaults.\n'); end
pop_editoptions( 'option_computeica', 1,'option_storedisk',0);

% commit ALLEEG variable to workspace
[data g] = hlp_splitstruct(g,{'ALLEEG'});
arg_toworkspace(data);
clear data;

% check the dataset(s)
ALLEEG = eeg_checkset(ALLEEG);

if procChans
    for cond=1:length(ALLEEG)
        % insert 'fake' ICA solution (Identity matrices) and copy EEG data into icaact
        [ALLEEG(cond).icaweights ALLEEG(cond).icasphere ALLEEG(cond).icawinv] = deal(eye(size(ALLEEG(cond).data,1)));
        ALLEEG(cond).icaact = ALLEEG(cond).data;
    end
end

for cond=1:length(ALLEEG)
    if ~isfield(ALLEEG(cond),'CAT') || ~isfield(ALLEEG(cond).CAT,'curComps')
        ALLEEG(cond).CAT.curComps = 1:size(ALLEEG(cond).icaweights,1);
    end
    if  ~isfield(ALLEEG(cond),'CAT') || ~isfield(ALLEEG(cond).CAT,'MODEL')
        ALLEEG(cond).CAT.MODEL = [];
    end
end


for cond=1:length(ALLEEG)
    
    g.EEG = ALLEEG(cond);
    
    if g.verb, fprintf('\nPre-processing condition %s\n', g.EEG.condition); end
    
    
    if isempty(g.EEG.icaact)
        g.EEG = hlp_icaact(g.EEG);
    end
    
    if g.backupOriginalData
        EEGbackup = g.EEG;  % save temporary copy of original dataset (TODO: clumsy, fix this)
    end
    
    % overwrite the data with ica activations
    % (this is a hack to allow us to use EEGLAB data selection routines)
    %     g.EEG.data = g.EEG.icaact;
    
    % preprocess the data
    % -------------------------------
    g  = hlp_preprocess(g);
    
    % store copy of processed data
    if g.verb
        fprintf('Storing processed data in g.EEG.CAT.srcdata...\n'); end
    
    % select components
    g.EEG.CAT.srcdata = g.EEG.icaact;
    g.EEG = pre_selectcomps('EEG',g.EEG,g.selectComponents,'verb',g.verb);
    g.EEG.CAT.nbchan = size(g.EEG.CAT.srcdata,1);
    
    % now also process the original dataset without normalization
    if g.backupOriginalData
        cfg = g;
        cfg.EEG     = EEGbackup;
        cfg.normalize.arg_selection = false;
        output      = hlp_preprocess(cfg);
        g.EEG.data  = output.EEG.data;
    end
    
    % estimate amplitude envelopes of selected data within a
    % specified frequency band
    if g.aamp.arg_selection
        g.EEG = est_aamp('EEG',g.EEG,g.aamp,'verb',g.verb);
    end
    
    % preserve times in case eeg_checkset deletes them (when num_trials=1)
    times = g.EEG.times;
    
    g.EEG.icaact        = [];  % force recompute of icaact
    EEGprep(cond)       = eeg_checkset(g.EEG);
    EEGprep(cond).times = times;
    args(cond)          = rmfield(g,'EEG');
    
    EEGprep(cond).CAT.configs.prepData = args(cond);
    EEGprep(cond).CAT.configs.prepData.arg_direct = 1;
    
    % force double precision
    if ~strcmpi(class(EEGprep(cond).CAT.srcdata),'double')
        EEGprep(cond).CAT.srcdata = double(EEGprep(cond).CAT.srcdata);
    end
    
end % next condition







function g = hlp_preprocess(g)

if ischar(g.normalize)
    g.normalize = {g.normalize};
end



% overwrite the data with ica activations
% (this is a hack to allow us to use EEGLAB data selection
% routines)
%         databackup = g.EEG.data;
%         g.EEG.data = g.EEG.icaact;


% select trials
if ~isempty(g.newtrials)
    if g.verb, disp('Selecting trials...'); end
    g.EEG = pop_select2(g.EEG,'trial',g.newtrials,'sorttrial','off');
    %         g.EEG.CAT.pre.newtrials = g.newtrials;
    if g.verb, fprintf('Done!\n'); end
end

% re-epoch data
if ~isempty(g.newtlims)
    if g.verb, disp(['Updating time limits to ' num2str(g.newtlims)]); end
    g.EEG = pop_select2(g.EEG,'time',g.newtlims);
    %         [dummy g.eventp] = min(abs(g.EEG.times));
    %         if g.EEG.pnts*g.EEG.srate > g.newtlims(end)
    %             fprintf('WARNING! endp=%d exceeds total number of points (%d). Updating endp to %d\n', ...
    %                 g.EEG.pnts, g.EEG.pnts, g.EEG.pnts);
    %
    %         end
    %         if g.varg.winlen > g.EEG.pnts/g.EEG.srate
    %             fprintf('WARNING! winlen=%0.1f s exceeds new epoch length (%0.1f s). Updating winlen to %0.1f s\n', ...
    %                 g.varg.winlen, g.EEG.pnts/g.EEG.srate, g.EEG.pnts/g.EEG.srate);
    %             g.varg.winlen = g.EEG.pnts/g.EEG.srate;
    %         end
    %         g.varg.endp = g.EEG.pnts;
    if g.verb, fprintf('Done!\n'); end
end

% detrend or center data
if g.detrend.arg_selection
    g.EEG = pre_detrend('EEG',g.EEG,g.detrend,'verb',g.verb);
end

%         filter data
if ~isempty(g.filter)
    if g.verb, fprintf('Filtering [%0.1f-%0.1f] Hz\n', g.filter(1), g.filter(2)); end
    
    if g.filter(2), g.EEG = pop_eegfilt(g.EEG,0,g.filter(2)); end % lowpass
    if g.filter(1), g.EEG = pop_eegfilt(g.EEG,g.filter(1),0); end % highpass
    
    %         g.EEG.CAT.pre.filtered = g.filter;
    if g.verb, fprintf('Done!\n'); end
end


% resample
if g.resample
    if g.verb, disp(['Downsampling to ' num2str(g.resample) ' Hz']); end
    
    g.EEG = pop_resample(g.EEG,g.resample);
    %         g.EEG.CAT.pre.resampled = srorig;
    if g.verb, fprintf('Done!\n'); end
end


% difference
if g.diff.arg_selection
    g.EEG = pre_diffData('EEG',g.EEG,g.diff,'verb',g.verb);
end

% remove bad segments of data
if isfield(g,'badsegments') && ~isempty(g.badsegments)
    for seg=1:size(g.badsegments,1)
        if g.verb, fprintf('Setting interval [%1.2f %1.2f] to NaN\n',g.badsegments(1),g.badsegments(2)); end
        [dummy pnts(1)] = min(abs(g.EEG.times-g.badsegments(1)*1000));
        [dummy pnts(2)] = min(abs(g.EEG.times-g.badsegments(2)*1000));
        g.EEG.data(:,pnts(1):pnts(2),:) = NaN;
    end
    if g.verb, fprintf('Done!\n'); end
end

g.EEG = eeg_checkset(g.EEG,'ica');

% normalize ica activations
if g.normalize.arg_selection
    g.EEG.icaact = pre_normData('data',g.EEG.icaact,g.normalize,'verb',g.verb);
end
