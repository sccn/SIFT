
function [MODEL cfg] = est_fitMVAR(varargin)
%
% Fit (adaptive) multivariate autoregressive model to EEG data. See [1] for
% details on VAR model fitting and implementations.
%
%
% Output:
%
%   MODEL structure with
%       .MODEL          (numvars x coeffs) matrix of VAR coefficients
%       .PE             (numvars x coeffs) prediction error (noise covariance) coefficients
%       .algorithm      string denoting algorithm used for estimation
%       .modelclass     string denoting model class (here, 'mvar')
%
% See Also: pop_est_fitMVAR(), pop_pre_prepData()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapters 3,6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift/
%
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

verb = arg_extract(varargin,{'verb','VerbosityLevel'},[],2);

g = arg_define([0 1],varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'EEGLAB dataset'), ...
    arg_subswitch({'algorithm','Algorithm'},hlp_getMVARalgorithms('defaultNameOnly'), ...
    hlp_getMVARalgorithms, ...
    {'Select a model fitting algorithm.',hlp_buildMVARHelpText}, ...
    'cat','Modeling Parameters','suppress',{'ModelOrder','OrderSelector','Verbosity','InitialState'}), ...
    arg({'morder','ModelOrder','modelOrder'},10,[1 Inf],'VAR model order.','cat','Modeling Parameters'), ...
    arg_nogui({'winStartIdx','WindowStartIndices'},[],[],'Starting indices for windows. This is a vector of sample points (start of windows) at which to estimate windowed VAR model','cat','Modeling Parameters'), ...
    arg({'winlen','WindowLength'},0.5,[eps Inf],'Sliding window length (sec)','cat','Modeling Parameters'), ...
    arg({'winstep','WindowStepSize'},0.03,[eps Inf],'Window step size (sec)','cat','Modeling Parameters'), ...
    arg({'epochTimeLims','EpochTimeLimits'},[],[],'Sub-epoch time limits (sec). This is relative to event time (e.g. [-1 2]). Default is the full epoch time range','cat','Modeling Parameters'), ...
    arg({'prctWinToSample','WindowSamplePercent'},100,[1 100],'Percent of windows to sample','cat','Modeling Parameters'), ...
    arg_subtoggle({'normalize','NormalizeData'},[],@pre_normData,'Z-normalize data within windows. Note this is not recommended for short windows','cat','Window Preprocessing','suppress',{'verb'}), ...
    arg_subtoggle({'detrend','Detrend'},{}, ...
    {arg({'method','DetrendingMethod'},'constant',{'linear','constant'},{'Detrend data within each window.', ...
    sprintf(['\n' ...
    'Linear: removes the least-squares fit of a straight line.\n' ...
    'Constant: removes the mean from each trial (centering)' ...
    ]) ...
    } ...
    )},'Detrend or center each time window','cat','Window Preprocessing'), ...
    arg({'timer','Timer'},false,[],'Activate timer. Times are stored in EEG.CAT.MODEL.timeelapsed'), ...
    arg({'verb','VerbosityLevel'},verb,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );

% commit EEG variable to workspace
[data g] = hlp_splitstruct(g,{'EEG'});
arg_toworkspace(data);
clear data;

% do some error-checking
if isempty(g.epochTimeLims)
    g.epochTimeLims = [EEG.xmin EEG.xmax]; end
if ~(all(g.epochTimeLims>=EEG.xmin) && all(g.epochTimeLims<=EEG.xmax))
    error('Epoch time limits must be within the range [%0.3g %0.3g]',EEG.xmin,EEG.xmax); end
if isempty(g.morder) || length(g.morder)>1
    error('invalid entry for field ''morder.'' Make sure the model order is a single integer.'); end

if nargout > 1, cfg = g; end

%% do some adjustments to parameters

% ensure we are using the right model order
g.algorithm.morder = g.morder;

if rem(g.winstep,1/EEG.srate)
    if g.winstep<1/EEG.srate
        g.winstep = 1/EEG.srate;
    else
        % adjust step size to nearest multiple of sampling interval
        g.winstep = g.winstep-rem(g.winstep,1/EEG.srate);
    end
    
    if g.verb,
        fprintf('Adjusting window step size to nearest multiple of sampling interval\n');
        fprintf('\tstep size is now %0.5g sec\n',g.winstep);
    end
end

if rem(g.winlen,1/EEG.srate)
    if g.winlen<1/EEG.srate
        g.winlen = 1/EEG.srate;
    else
        % adjust window length to nearest multiple of sampling interval
        g.winlen = g.winlen-rem(g.winlen,1/EEG.srate);
    end
    
    if g.verb,
        fprintf('Adjusting window length to nearest multiple of sampling interval\n');
        fprintf('\twindow length is now %0.5g sec\n',g.winlen);
    end
end

tidx = getindex(EEG.CAT.times,g.epochTimeLims*1000);
if ~all(isequal(EEG.CAT.times(tidx),g.epochTimeLims*1000))
    
    g.epochTimeLims = EEG.CAT.times(tidx)/1000;
    
    if g.verb
        fprintf('Adjusting epoch time limits to match sampling interval\n');
        fprintf('\tepoch limits are now [%0.5g, %0.5g] sec\n',g.epochTimeLims(1),g.epochTimeLims(2));
    end
end

winLenPnts  = round(g.winlen*EEG.srate); % window size in points
winStepPnts = round(g.winstep*EEG.srate);

if isempty(g.winStartIdx)
    % starting point of each window (points)
    g.winStartIdx  = tidx(1):winStepPnts:(tidx(2)-winLenPnts)+1;
    %g.winStartIdx  =  round((double(g.epochTimeLims(1):g.winstep:g.epochTimeLims(2)-g.winlen)*EEG.srate)+1;
end

if g.prctWinToSample<100
    % randomly select percentage of windows to work with
    randwin = randperm(length(g.winStartIdx));
    randwin = sort(randwin(1:ceil(length(g.winStartIdx)*g.prctWinToSample/100)));
    g.winStartIdx = g.winStartIdx(randwin);
end

numWins   = length(g.winStartIdx);

%% Prepare data for model fitting

% initialize results arrays
[AR PE RC mu th]  = deal(cell(1,numWins));

if g.verb
    waitbarstring = sprintf(['Fitting VAR[%d] model using %s algorithm\n' ...
        fastif(isempty(EEG.condition),'',['Condition: ' EEG.condition])], ...
        g.morder, num2str(g.algorithm.arg_selection));
    h=waitbar(0,waitbarstring);
end

if g.normalize.arg_selection
    % normalize each window separately
    if g.verb, fprintf('I will normalize each window across %s...\n',g.normalize.method{:}); end
%     for t=1:numWins
%         winpnts = g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1;
%         EEG.CAT.srcdata(:,winpnts,:) = pre_normData(EEG.CAT.srcdata(:,winpnts,:),'Method',g.normalize.method,'verb',0);
%     end
end

if g.detrend.arg_selection
    % detrend each window separately
    if g.verb, fprintf('%s detrending each window...\n',firstcaps(g.detrend.method)); end
    for ch=1:EEG.CAT.nbchan
        EEG.CAT.srcdata(ch,:,:) = locdetrend_siftmod(squeeze(EEG.CAT.srcdata(ch,:,:)), ...
            EEG.srate,[g.winlen g.winstep],g.detrend.method);
    end
end

if g.timer
    timeElapsed = nan(1,numWins);
else
    timeElapsed = [];
end

%% Main loop: fit MVAR models to each window
for t=1:numWins
    
    if g.timer, tic; end
    
    % get indices of all data points (samples) in current window
    winpnts = g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1;
    
    % select the data chunk
    datachunk = EEG.CAT.srcdata(:,winpnts,:);
    
    if g.normalize.arg_selection
        % normalize the data chunk
        datachunk = pre_normData(datachunk,'Method',g.normalize.method,'verb',0);
    end

    % execute the model-fitting algorithm
    algFcnName = hlp_getMVARalgorithms('mfileNameOnly',g.algorithm.arg_selection);
    switch nargout(algFcnName)
        case 2
            [AR{t} PE{t}] = feval(algFcnName, ...
                                  'data',datachunk, ...
                                   g.algorithm);
        case 3
            [AR{t} PE{t} argsout] = feval(algFcnName, ...
                                   'data',datachunk,...
                                    g.algorithm);
            if isstruct(argsout)
                % store contents of argsout fields in cell array at index t
                % e.g. fieldname{t} = argsout.(fieldname)
                fnames = fieldnames(argsout);
                for k=1:length(fnames)
                    eval([fnames{k} '{t}=argsout.(''' fnames{k} ''');']);
                end
            end
        otherwise
            error('SIFT:est_fitMVAR:badAlgArgs','%s must output either 2 or 3 arguments',algFcnName);
    end
        
        
    if g.verb
        try
            waitbar(t/numWins,h,...
                sprintf(' %s (%d/%d)',waitbarstring,t,numWins));
        catch err
            if strcmp(err.identifier,'MATLAB:waitbar:InvalidArguments');
                MODEL = [];
                return;
            end
        end
    end
    
    if g.timer, timeElapsed(t) = toc; end
end

%% Do some cleanup
clear(algFcnName);
clear('-regexp','mvar_*');
if g.verb, close(h); end

%% Construct MODEL object
switch lower(g.algorithm.arg_selection)
    case 'group lasso dal/scsa'
        %     MODEL.ww = ww;
        MODEL.lambda = g.algorithm.dal_args.lambda;
    case 'group lasso admm'
        MODEL.lambda = g.algorithm.admm_args.lambda;
        MODEL.rho    = g.algorithm.admm_args.rho;
        MODEL.Alpha  = g.algorithm.admm_args.Alpha;
end

MODEL.AR = AR;
MODEL.PE = PE;
MODEL.RC = RC;
MODEL.mu = mu;
MODEL.th = th;
MODEL.winStartTimes = (g.winStartIdx-1)/EEG.srate; %EEG.CAT.times(g.winStartIdx)/1000;
MODEL.morder        = g.morder;
MODEL.winstep       = g.winstep;
MODEL.winlen        = g.winlen;
MODEL.algorithm     = g.algorithm.arg_selection;
MODEL.modelclass    = 'mvar';
MODEL.timeelapsed   = timeElapsed;
MODEL.normalize     = g.normalize;
MODEL.modelapproach = 'Segmentation VAR';





function algs = hlp_getMVARalgorithms(varargin)
% return a cell array of MVAR algorithm names and associated function
% handles. The format is
% algs = {{name1 func_handle1} {name2 func_handle2} ... }
% if defaultNameOnly = true, then we only return name1 as a string
%
% hlp_getMVARalgorithms('defNameOnly') returns only the name of the 
%  default algorithm (as a string)
%
% hlp_getMVARalgorithms('mfileNameOnly',algName) where algName is a string
%   with the human-readable name of a valid algorithm, returns the m-file 
%   name of the function implementing the specified approach
% 
% NOTES:
%
% If the preamble text of the function (c.f. hlp_getFcnPreambleText())
% begins with the line
%
% Algorithm: <algorithm_name>
%
% then <algorithm_name> is returned as the human-readable algorithm name. 
% Otherwise, if this line cannot be found, the function name is returned.

persistent mvarAlgs;

defaultNameOnly = false;
algName = '';

if nargin==1 && strcmp(varargin{1},'defaultNameOnly')
    defaultNameOnly = true;
end

if nargin==2 && ismember('mfileNameOnly',varargin);
    algName = varargin{2};
    if ~ischar(algName)
        error('SIFT:est_fitMVAR:hlp_getMVARalgorithms:badInput', ...
              'Bad argument pair for ''mfileNameOnly'''); 
    end
end

if ~isempty(mvarAlgs)
    % only do the check once (it's time-consuming)
    % ... retrive "cached" version
    algs = mvarAlgs;
else
    
    algs = {};
        
    % get the names of all mvar_* functions in the /est/mvar folder
    siftroot = fileparts(which('StartSIFT'));
    fpath    = [siftroot filesep 'est' filesep 'mvar' filesep];
    mvarFcns = wildcardsearch(fpath,'*mvar_*.m',true,true);
    mvarFcns = regexprep(mvarFcns,['.*mvar' filesep],'');
    mvarFcns = regexprep(mvarFcns,'\.m','');
    
    % cycle through the list of algorithms
    for k=1:length(mvarFcns)
        
        % initialize variable which will determine whether 
        % we exclude this algorithm (for instance if one or
        % more dependencies are missing)
        skipAlgorithm = false; 
        
        % get the help text (H1) for the algorithm entry function
        preText = hlp_getFcnPreambleText(mvarFcns{k});
        
        % extract the human-readable algorithm name from the text 
        % following the 'Algorithm:' header
        preText = strtrim(preText);
        algName = regexpi(preText,'Algorithm\s*[:]?\s*([^\n]*)','tokens');
        if ~isempty(algName)
            algName = strtrim(algName{1}{1});
        else
            algName = mvarFcns{k};
        end
        
        % search if there are any key dependencies (functions) missing from
        % the path. Key dependencies (m-file names) are listed in the 
        % 'Dependencies:' block of the function HelpText preamble.
        deplist = regexpi(preText,'Dependencies\s*[:]?\s*([^\n]*)','tokens');
        if ~isempty(deplist)
            % parse the list of function dependencies 
            % and extract m-file/function names
            deplist = regexp(deplist{1}{1},'\S*[^ \(\),\s]','match');
            deplist = strtrim(deplist); % just in case

            for dep=1:length(deplist)
                % check if the dependency function exists
                if ~exist(deplist{dep},'file')
                    % ... a critical dependency does not exist
                    % so we have to exclude this algorithm
                    skipAlgorithm = true;
                    break;
                end
            end
        end
        
        if ~skipAlgorithm
            algs{end+1} = {algName str2func(mvarFcns{k})};
        end
    end
    
    mvarAlgs = algs;
end

% return only the name of the first available algorithm
if defaultNameOnly
    algs = 'Vieira-Morf';
%     algs = algs{1}{1};
elseif ~isempty(algName)
    % get the function handle matching the desired algorithm name
    tmp = cellfun(@(x) fastif(isequal(x{1},algName),x{2},''),algs,'UniformOutput',false);
    tmp(cellfun(@isempty,tmp))=[];
    if isempty(tmp)
        error('SIFT:est_fitMVAR:hlp_getMVARalgorithms:badAlgorithmName', ...
              'Unknown algorithm ''%s''',algName);
    else
        % return function name as a string
        algs = char(tmp{1});
    end
end



function hlpText = hlp_buildMVARHelpText()
% return a formatted help text string which contains help text definitions
% for each of the var modeling functions

% first get the list of function names
algs = hlp_getMVARalgorithms;

hlpText = '';

% for each algorithm...
for k=1:length(algs)
    % ... get the human-readable name for the algorithm
    algFullName = algs{k}{1};
    % ... get the function m-file name ...
    algFcnName = hlp_getMVARalgorithms('mfileNameOnly',algFullName);
    % ... get the preamble for the function ...
    preamble = hlp_getFcnPreambleText(algFcnName);
    % ... and insert into the hlpText
    hlpText = sprintf([hlpText '\n\n' ...
                       '---------------------------------------\n' ...
                       algFullName '\n' ...
                       '---------------------------------------\n' ...
                       preamble
                       ]);
end