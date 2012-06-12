
function [MODEL params] = est_fitMVAR(varargin)
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
    arg_subswitch({'algorithm','Algorithm'},hlp_getMVARalgorithms(true), ...
    hlp_getMVARalgorithms(false), ...
    {'Select a model fitting algorithm.', ...
        sprintf([...
            '\n' ...
            '\n' ...
            '------------------\n' ...
            'Vieira-Morf:\n' ...
            '------------------\n' ...
            'Unconstrained VAR modeling via Vieira-Morf Maximum Entropy algorithm.\n', ...
            '\n' ...
            'References and code:\n', ...
            '[1] A. Schlogl, Comparison of Multivariate Autoregressive Estimators. Signal processing, Elsevier B.V.\n', ...
            '[2] S.L. Marple "Digital Spectral Analysis with Applications" Prentice Hall, 1987.\n', ...
            '\n' ...
            '------------------\n' ...
            'GroupLasso_ADMM:\n' ...
            '------------------\n' ...
            'Sparse VAR modeling via Group Lasso.\n', ...
            'This option estimates sparse VAR coefficients using the Alternating Direction Method of Multipliers (ADMM).\n', ...
            '\n' ...
            'References and code:\n', ...
            '[1] Boyd, Parikh, Chu, Pelato and Eckstein et al. Foundations and Trends in Machine Learning 3(1):1-122,2011.\n', ...
            '[2] http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html' ...
            '\n' ...
            ])},'cat','Modeling Parameters','suppress',{'ModelOrder','OrderSelector','Verbosity','InitialState'}), ...
    arg({'morder','ModelOrder','modelOrder'},10,[1 Inf],'VAR model order.','cat','Modeling Parameters'), ...
    arg_nogui({'winStartIdx','WindowStartIndices'},[],[],'Starting indices for windows. This is a vector of sample points (start of windows) at which to estimate windowed VAR model','cat','Modeling Parameters'), ...
    arg({'winlen','WindowLength'},0.5,[eps Inf],'Sliding window length (sec)','cat','Modeling Parameters'), ...
    arg({'winstep','WindowStepSize'},0.03,[eps Inf],'Window step size (sec)','cat','Modeling Parameters'), ...
    arg({'epochTimeLims','EpochTimeLimits'},[],[],'Sub-epoch time limits (sec). This is relative to event time (e.g. [-1 2]). Default is the full epoch time range','cat','Modeling Parameters'), ...
    arg({'prctWinToSample','WindowSamplePercent'},100,[1 100],'Percent of windows to sample','cat','Modeling Parameters'), ...
    arg_subtoggle({'normalize','NormalizeData'},[],@pre_normData,'Z-normalize data within windows. Note this is not recommended for short windows','cat','Window Preprocessing'), ...
    arg_subtoggle({'detrend','Detrend'},{}, ...
    {arg({'method','DetrendingMethod'},'constant',{'linear','constant'},{'Detrending options.', ...
    sprintf(['\n' ...
    'Linear: removes the least-squares fit of a straight line from each trial.\n' ...
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

% % initialize defaults for scsa
% scsadef = struct('lambda',100,'shrink_diagonal',true,'initAR',[],'loss','hs');
% if ~isfield(g,'scsa'), g.scsa = scsadef; end
% fnames = fieldnames(g.scsa);
% for f=1:length(fnames)
%     scsadef.(fnames{f}) = g.scsa.(fnames{f});
% end
% g.scsa = scsadef;
%
% % initialize defaults for glADMM
% gladmmdef = struct('lambda',[],'initAR',[],'rho',1.0,'alpha',1.0,'verb',0,'max_iter',1000);
% if ~isfield(g,'gladmm'), g.gladmm = gladmmdef; end
% fnames = fieldnames(g.gladmm);
% for f=1:length(fnames)
%     gladmmdef.(fnames{f}) = g.gladmm.(fnames{f});
% end
% g.gladmm = gladmmdef;
%
% % initialize defaults for ridge regression
% ridgedef = struct('lambda',0.1,'scaled',1);
% if ~isfield(g,'ridge'), g.ridge = ridgedef; end
% fnames = fieldnames(g.ridge);
% for f=1:length(fnames)
%     ridgedef.(fnames{f}) = g.ridge.(fnames{f});
% end
% g.ridge = ridgedef;


if nargout > 1, params = g; end

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

tidx = getindex(EEG.times,g.epochTimeLims*1000);
if ~all(isequal(EEG.times(tidx),g.epochTimeLims*1000))
    
    g.epochTimeLims = EEG.times(tidx)/1000;
    
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
    if g.verb, fprintf('Normalizing each window across %s...\n',g.normalize.method{:}); end
    for t=1:numWins
        winpnts = g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1;
        EEG.CAT.srcdata(:,winpnts,:) = pre_normData(EEG.CAT.srcdata(:,winpnts,:),'Method',g.normalize.method,'verb',0);
    end
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

% % check if algorithms exist
% exist_mvarpll   = logical(exist('mvar_pll','file'));
% exist_armorf    = logical(exist('armorf','file'));
% exist_arfit     = logical(exist('arfit','file'));
% exist_scsa      = logical(exist('dalhsgl','file'));


%% Main loop: fit MVAR models to each window
for t=1:numWins
    
    if g.timer, tic; end
    
    % get indices of all data points (samples) in current window
    winpnts = g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1;
    
    % select the appropriate algorithm and execute
    switch lower(g.algorithm.arg_selection)
        case 'vieira-morf-pll'
            if 0 %exist_mvarpll
                
                if t==1
                    % permute to [time x chans x trials]
                    EEG.CAT.srcdata = permute(EEG.CAT.srcdata,[2 1 3]);
                end
                
                data = squeeze(EEG.CAT.srcdata(winpnts,:,:));
                % make 2-D with ModelOrder+2 NaNs between trials
                data = nanpad(data,g.morder);
                [AR{t},RC{t},PE{t}] = mvar_pll(gdouble(data), g.morder, 2);
            else
                error('mvar_pll.m not found! mvar_pll option not available!');
            end
        case 'vieira-morf-ding'
            if 0 %exist_armorf
                [AR{t},PE{t}] = armorf(reshape(EEG.CAT.srcdata(:,winpnts,:), ...
                    [EEG.CAT.nbchan,EEG.trials*winLenPnts]),...
                    EEG.trials,winLenPnts,g.morder);
                AR{t} = -AR{t}; % convert to default mvar format
            else
                error('armorf.m not found! vieira-morf-ding option unavailable');
            end
        case 'vieira-morf'
            
            if t==1
                % permute to [time x chans x trials]
                EEG.CAT.srcdata = permute(EEG.CAT.srcdata,[2 1 3]);
            end
            
            data = squeeze(EEG.CAT.srcdata(winpnts,:,:));
            % make 2-D with ModelOrder+2 NaNs between trials
            data = nanpad(data,g.morder);
            [AR{t},RC{t},PE{t}] = mvar_vieiramorf('data',data, g.algorithm);
        case 'arfit'
            RC{t} = [];
            [mu{t}, AR{t}, PE{t} sbc, fpe, th{t}] = mvar_arfit('data',permute(EEG.CAT.srcdata(:,winpnts,:),[2 1 3]),g.algorithm);
            
        case 'group lasso dal/scsa'
            if t==1
                % initialization block
                %                     if iscell(g.algorithm.AR0)
                %                         AR0 = g.algorithm.AR0{t};
                %                     else
                %                         AR0 = g.algorithm.AR0;
                %                     end
                
                if ~strcmpi(class(EEG.CAT.srcdata),'double')
                    % convert data to double precision (necessary for kron)
                    EEG.CAT.srcdata = double(EEG.CAT.srcdata);
                end
            end
            
            [AR{t} PE{t}] = mvar_dalSCSA('data',EEG.CAT.srcdata(:,winpnts,:),g.algorithm);
            RC{t} = [];
        case 'ridge regression'
            [AR{t} PE{t}] = mvar_ridge('data',EEG.CAT.srcdata(:,winpnts,:),g.algorithm);
            RC{t} = [];
        case 'group lasso admm'
            if t==1
                % intialization block
                %                 if iscell(g.algorithm.AR0)
                %                     AR0 = g.algorithm.AR0{t};
                %                 else
                %                     AR0 = g.algorithm.AR0;
                %                 end
                if ~strcmpi(class(EEG.CAT.srcdata),'double')
                    % convert data to double precision (necessary for kron)
                    EEG.CAT.srcdata = double(EEG.CAT.srcdata);
                end
            end
            
            [AR{t} PE{t}] = mvar_glADMM('data',EEG.CAT.srcdata(:,winpnts,:),g.algorithm);
            RC{t} = [];
        otherwise
            if 0 %exist('mvar','file')
                
                if t==1
                    % permute to [time x chans x trials]
                    EEG.CAT.srcdata = permute(EEG.CAT.srcdata,[2 1 3]);
                end
                % one of the other mvar modes...
                data = squeeze(EEG.CAT.srcdata(winpnts,:,:));
                % make 2-D with ModelOrder+2 NaNs between trials
                data = nanpad(data,g.morder);
                [AR{t}, RC{t}, PE{t}] = mvar(data,g.morder,str2double(g.algorithm.arg_selection));
            else
                help 'est_fitMVAR';
                error('unknown algorithm (%s)',g.algorithm.arg_selection);
            end
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

switch lower(g.algorithm.arg_selection)
    case 'scsa'
        % clear persistent variables
        clear mvar_dalSCSA;
        %     MODEL.ww = ww;
        MODEL.lambda = g.algorithm.dal_args.lambda;
    case 'group lasso admm'
        clear mvar_glADMM
        
        MODEL.lambda = g.algorithm.admm_args.lambda;
        MODEL.rho    = g.algorithm.admm_args.rho;
        MODEL.Alpha  = g.algorithm.admm_args.Alpha;
end

if g.verb, close(h); end


% construct MODEL object
MODEL.AR = AR;
MODEL.PE = PE;
MODEL.RC = RC;
MODEL.mu = mu;
if strcmpi(g.algorithm.arg_selection,'arfit')
    MODEL.th = th;
end
MODEL.winStartTimes = (g.winStartIdx-1)/EEG.srate; %EEG.times(g.winStartIdx)/1000;
MODEL.morder        = g.morder;
MODEL.winstep       = g.winstep;
MODEL.winlen        = g.winlen;
MODEL.algorithm     = g.algorithm.arg_selection;
MODEL.modelclass    = 'mvar';
MODEL.timeelapsed   = timeElapsed;
MODEL.normalize     = g.normalize;




function algs = hlp_getMVARalgorithms(defaultNameOnly)
% return a cell array of MVAR algorithm names and associated function
% handles. The format is
% algs = {{name1 func_handle1} {name2 func_handle2} ... }
% if defaultNameOnly = true, then we only return name1 as a string
persistent mvarAlgs;

if isempty(mvarAlgs)
    % only do the check once (it's time-consuming)
    
    algs = {};
    
    if exist('mvar_vieiramorf','file')
        algs{end+1} = {'Vieira-Morf'           @mvar_vieiramorf};   end
    if exist('arfit','file')
        algs{end+1} = {'ARfit'                 @mvar_arfit};        end
    if exist('dalhsgl','file')
        algs{end+1} =  {'Group Lasso DAL/SCSA'  @mvar_dalSCSA};      end
    if exist('mvar_glADMM','file')
        algs{end+1} =  {'Group Lasso ADMM'      @mvar_glADMM};       end
    if exist('ridge','file')
        algs{end+1} =  {'Ridge Regression'      @mvar_ridge};        end
    
    mvarAlgs = algs;
else
    
    % retrive "cached" version
    algs = mvarAlgs;
end

if defaultNameOnly
    % return only the name of the first available algorithm
    algs = algs{1}{1};
end
