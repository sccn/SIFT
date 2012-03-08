
function [MODEL params] = est_fitMVAR(EEG,typeproc,varargin)
%
% Fit (adaptive) multivariate autoregressive model to EEG data. See [1] for
% details on VAR model fitting and implementations.
%
% Input:
%
%   EEG                Preprocessed EEG structure.
%   typeproc           Reserved for future use. Use 0
%
% Optional:
%
% 'algorithm',         string denoting which algorithm to use for model
%                       fitting ('vieira-morf','arfit','scsa')
% 'winStartIdx'        vector of sample points (start of windows) at which to estimate windowed VAR model
% 'morder',            VAR model order
% 'winlen',            window length (sec)
% 'winstep',           window step size (sec)
% 'epochTimeLims',     time range to analyze (sec) where 0 = start of the epoch
% 'prctWinToSample',   percent of time windows to randomly select
% 'verb',              verbosity level (0=no output, 1=text, 2=gui)
% 'timer'              estimate time required to fit model
% 'normalize'          cell array containing one or more of
%                      {'temporal', 'ensemble'}. This performs ensemble
%                      normalization or temporal normalization (or both)
%                      within each window
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


if nargin<3
    help 'est_fitMVAR';
end

var = hlp_mergeVarargin(varargin{:});
g = finputcheck(var, hlp_getDefaultArglist('est'), 'est_fitMVAR','ignore','quiet');
if ischar(g), error(g); end
if isempty(g.epochTimeLims), g.epochTimeLims = [0 EEG.xmax]; end
if isempty(g.morder) || length(g.morder)>1, error('invalid entry for field ''morder.'' Make sure the model order is a single integer.'); end

% initialize defaults for scsa
scsadef = struct('lambda',100,'shrink_diagonal',true,'initAR',[],'loss','hs');
if ~isfield(g,'scsa'), g.scsa = scsadef; end
fnames = fieldnames(g.scsa);
for f=1:length(fnames)
    scsadef.(fnames{f}) = g.scsa.(fnames{f});
end
g.scsa = scsadef;

% initialize defaults for glADMM
gladmmdef = struct('lambda',[],'initAR',[],'rho',1.0,'alpha',1.0,'verb',0,'max_iter',1000);
if ~isfield(g,'gladmm'), g.gladmm = gladmmdef; end
fnames = fieldnames(g.gladmm);
for f=1:length(fnames)
    gladmmdef.(fnames{f}) = g.gladmm.(fnames{f});
end
g.gladmm = gladmmdef;

% initialize defaults for ridge regression
ridgedef = struct('lambda',0.1,'scaled',1);
if ~isfield(g,'ridge'), g.ridge = ridgedef; end
fnames = fieldnames(g.ridge);
for f=1:length(fnames)
    ridgedef.(fnames{f}) = g.ridge.(fnames{f});
end
g.ridge = ridgedef;


if nargout > 1, params = g; end

if isempty(g.algorithm)
    g.algorithm = 'vieira-morf';
end

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

winLenPnts = floor(g.winlen*EEG.srate); % window size in points
if isempty(g.winStartIdx)
    % starting point of each window (points)
    g.winStartIdx  = round((g.epochTimeLims(1):g.winstep:g.epochTimeLims(2)-g.winlen)*EEG.srate)+1;
end
if g.prctWinToSample<100
    % randomly select percentage of windows to work with
    randwin = randperm(length(g.winStartIdx));
    randwin = sort(randwin(1:ceil(length(g.winStartIdx)*g.prctWinToSample/100)));
    g.winStartIdx = g.winStartIdx(randwin);
end
numWins   = length(g.winStartIdx);

% initialize variables
% th contains confidence intervals from ARFIT (if used)
[AR PE RC th]  = deal(cell(1,numWins));

if g.verb
    h=waitbar(0,sprintf('fitting VAR[%d] model [mode=%s]\nCondition: %s ...', ...
        g.morder, num2str(g.algorithm),EEG.condition));
end

if ~isempty(g.normalize)
    % normalize each window separately
    if g.verb, fprintf('Normalizing each window across %s...\n',g.normalize{:}); end
    for t=1:numWins
        winpnts = g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1;
        EEG.CAT.srcdata(:,winpnts,:) = pre_normData(EEG.CAT.srcdata(:,winpnts,:),'Method',g.normalize);
    end
end

if ~isempty(g.detrend)
    % detrend each window separately
    if g.verb, fprintf('%s detrending each window...\n',firstcaps(g.detrend)); end
    for ch=1:EEG.CAT.nbchan
        EEG.CAT.srcdata(ch,:,:) = locdetrend_siftmod(squeeze(EEG.CAT.srcdata(ch,:,:)), ...
            EEG.srate,[g.winlen g.winstep],g.detrend);
    end
end

if g.timer
    timeElapsed = nan(1,numWins);
else
    timeElapsed = [];
end

if ~ismember(g.algorithm,{'arfit','vieira-morf-ding','scsa','ridge','gladmm'})
    EEG.CAT.srcdata = permute(EEG.CAT.srcdata,[2 1 3]);  % time x chans x trials
end

if isequal(g.algorithm,'scsa') && ~strcmpi(class(EEG.CAT.srcdata),'double')
    % convert data to double precision (necessary for kron)
    EEG.CAT.srcdata = double(EEG.CAT.srcdata);
end

exist_mvarpll   = logical(exist('mvar_pll','file'));
exist_armorf    = logical(exist('armorf','file'));
exist_arfit     = logical(exist('arfit','file'));
exist_scsa      = logical(exist('dalhsgl','file'));

for t=1:numWins
    
    if g.timer, tic; end
    % extract window and pad each trial with nans up to model order+2
    % Fit MODEL model up to g.morder
    winpnts = g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1;
    
    
    switch lower(g.algorithm)
        case 'vieira-morf-pll'
            if exist_mvarpll
                data = squeeze(EEG.CAT.srcdata(winpnts,:,:));
                data = nanpad(data,g.morder);
                [AR{t},RC{t},PE{t}] = mvar_pll(gdouble(data), g.morder, 2);
            else
                error('mvar_pll.m not found! mvar_pll option not available!');
            end
        case 'vieira-morf-ding'
            if exist_armorf
                [AR{t},PE{t}] = armorf(reshape(EEG.CAT.srcdata(:,winpnts,:), ...
                    [EEG.CAT.nbchan,EEG.trials*winLenPnts]),...
                    EEG.trials,winLenPnts,g.morder);
                AR{t} = -AR{t}; % convert to default mvar format
            else
                error('armorf.m not found! vieira-morf-ding option unavailable');
            end
        case 'vieira-morf'
            data = squeeze(EEG.CAT.srcdata(winpnts,:,:));
            data = nanpad(data,g.morder);
            [AR{t},RC{t},PE{t}] = mvar_vieiramorf(data, g.morder);
        case 'arfit'
            if exist_arfit
                RC{t} = [];
                [w, AR{t}, PE{t} sbc, fpe, th{t}] = arfit(permute(EEG.CAT.srcdata(:,winpnts,:),[2 1 3]), g.morder, g.morder,'zero');
            else
                error('arfit.m not found! Please download ARFIT. ARFIT option unavailable');
            end
        case 'scsa'
            if exist_scsa
                if iscell(g.scsa.initAR)
                    AR0 = g.scsa.initAR{t};
                else
                    AR0 = g.scsa.initAR;
                end
                [AR{t} PE{t}] = mvar_dalSCSA(EEG.CAT.srcdata(:,winpnts,:),g.morder,g.scsa.lambda,AR0,g.scsa.shrink_diagonal,g.scsa.loss,'solver','cg','maxiter',10,'tol',1e-3);
                RC{t} = [];
            else
                error('dalhsgl.m not found! Please download DAL 1.05 toolbox. SCSA option unavailable');
            end
        case 'ridge'
            [AR{t} PE{t}] = mvar_ridge(EEG.CAT.srcdata(:,winpnts,:),g.morder,g.ridge.lambda,g.ridge.scaled);
            RC{t} = [];
        case 'gladmm'
            if iscell(g.gladmm.initAR)
                AR0 = g.gladmm.initAR{t};
            else
                AR0 = g.gladmm.initAR;
            end
            
            [AR{t} PE{t}] = mvar_glADMM(EEG.CAT.srcdata(:,winpnts,:),g.morder,g.gladmm.lambda,g.gladmm.rho,g.gladmm.alpha,AR0,g.gladmm.verb,g.gladmm.max_iter);
            RC{t} = [];
        otherwise
            if exist('mvar','file')
                % one of the other mvar modes...
                data = squeeze(EEG.CAT.srcdata(winpnts,:,:));
                data = nanpad(data,g.morder);
                [AR{t}, RC{t}, PE{t}] = mvar(data,g.morder,str2double(g.algorithm));
            else
                help 'est_fitMVAR';
                error('unknown algorithm (%s)',g.algorithm);
            end
    end
    
    if g.verb
        try
            waitbar(t/numWins,h,...
                sprintf('fitting VAR[%d] model [mode=%s] (window %d/%d)\nCondition: %s...', ...
                g.morder,num2str(g.algorithm),t,numWins,EEG.condition));
        catch err
            if strcmp(err.identifier,'MATLAB:waitbar:InvalidArguments');
                MODEL = [];
                return;
            end
        end
    end
    
    if g.timer, timeElapsed(t) = toc; end
end

switch lower(g.algorithm)
    case 'scsa'
        % clear persistent variables
        clear mvar_dalSCSA
        %     MODEL.ww = ww;
        MODEL.lambda = g.scsa.lambda;
    case 'gladmm'
        clear mvar_glADMM
        
        MODEL.lambda = g.gladmm.lambda;
        MODEL.rho    = g.gladmm.rho;
        MODEL.alpha  = g.gladmm.alpha;
end

if g.verb, close(h); end

% construct MODEL object
MODEL.AR = AR;
MODEL.PE = PE;
MODEL.RC = RC;
if strcmpi(g.algorithm,'arfit')
    MODEL.th = th;
end
MODEL.winStartTimes = (g.winStartIdx-1)/EEG.srate;
MODEL.morder = g.morder;
MODEL.winstep = g.winstep;
MODEL.winlen = g.winlen;
MODEL.algorithm = g.algorithm;
MODEL.modelclass = 'mvar';
MODEL.timeelapsed = timeElapsed;
MODEL.normalize = g.normalize;
