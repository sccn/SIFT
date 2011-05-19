
function [MODEL params] = est_fitMVARKalman(EEG,typeproc,varargin)
%
% Fit multivariate autoregressive model to EEG data using a 
% Kalman Filter. This function is a wrapper for the mvaar() function from
% the TSA toolbox [2]. See [1] for additional details on VAR model fitting 
% and implementation.
%
% Input:
%
%   EEG                Preprocessed EEG structure.
%   typeproc           Reserved for future use. Use 0
%
% Optional:
%
%     'updatecoeff'        Kalman filter update coefficient {def: 0.001}
%     'updatemode'         Kalman filter noise update mode {def: 1}
%     'winStartIdx'        vector of sample points (start of windows) at which to estimate windowed VAR model
%     'morder'             VAR model order
%     'epochTimeLims'      time range to analyze (sec) where 0 = start of the epoch
%     'verb'               verbosity level (0=no output, 1=text, 2=gui)
%     'timer'              estimate time required to fit model
%     'downsampleFactor'   store VAR coefficient matrices only for every Nth
%                          timepoint, where N=downsampleFactor.
%
% Output:
%
%   MODEL structure with
%       .AR             (numvars x coeffs) matrix of VAR coefficients
%       .PE             (numvars x coeffs) prediction error (noise covariance) coefficients
%       .algorithm      string denoting algorithm used for estimation
%       .modelclass     string denoting model class (here, 'mvar')
%
% See Also: est_fitMVAR(), mvaar()
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual. Chapters 3,6. 
%     Available at: http://www.sccn.ucsd.edu/wiki/Sift
% [2] A. Schloegl. The time series analysis toolbox for octave and matlab, 
%     available online at http://www.dpmi.tu-graz.ac. at/~schloegl/matlab/tsa/ 
%     and http://cvs.sourceforge.net/ viewcvs.py/octave/octave-forge/extra/tsa/.
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
    help 'est_fitMVARKalman';
end


var = hlp_mergeVarargin(varargin{:});
g = finputcheck(var, hlp_getDefaultArglist('est'), 'est_fitMVAR','ignore');
if ischar(g), error(g); end
if isempty(g.epochTimeLims), g.epochTimeLims = [0 EEG.pnts/EEG.srate]; end
if isempty(g.morder) || length(g.morder)>2, error('invalid entry for field ''morder'''); end
 % combine structs, overwriting duplicates of g with 
if ~isfield(g,'updatecoeff'), g.updatecoeff = 0.001; end
if ~isfield(g,'updatemode'), g.updatemode = 1; end
if ~isfield(g,'downsampleFactor'), g.downsampleFactor = []; end;
if ~isfield(g,'constraints'), g.constraints = struct([]); end    % constraints.D and constraints.d are constraints of form Dx=d

if isempty(g.downsampleFactor)
    g.downsampleFactor = round(g.winstep*EEG.srate);
end

g.winstep = g.downsampleFactor/EEG.srate;


%     g = catstruct(g,gvar); clear g2;

if nargout > 1, params = g; end

g.winStartIdx  = 1:g.downsampleFactor:size(EEG.CAT.srcdata,2);
numWins   = length(g.winStartIdx);

[AR PE RC]  = deal(cell(1,numWins));

if g.timer
    timeElapsed = nan(1,numWins);
else
    timeElapsed = [];
end
if g.timer, tic; end

if size(EEG.CAT.srcdata,3)>1
    error('Kalman filtering cannot be used with multi-trial data');
end

MemoryLengthSamples = -1/log(1-g.updatecoeff);
TimeConstant = MemoryLengthSamples/EEG.srate;
if g.verb, 
    fprintf('The effective window length is approximately %0.5f seconds (%d samples)\n',TimeConstant,MemoryLengthSamples); 
    fprintf('Your step size is %0.3f ms\n',g.winstep*1000);
    if ~isempty(g.constraints)
        fprintf('Using constraints\n');
    end
end
    

% Fit MODEL model up to g.morder using Kalman filter
 
[nchs npnts] = size(EEG.CAT.srcdata);

[VAR,residuals,Kalman,C] = mvaar(permute(EEG.CAT.srcdata,[2 1]),g.morder,g.updatecoeff,g.updatemode,[],g.verb,g.downsampleFactor,g.constraints);

time = toc;

for t=2:size(VAR,1)
    % store VAR coefficients and noise covariance matrices
    AR{t-1} = reshape(VAR(t,:),nchs*g.morder,nchs).';
    PE{t-1}  = C(:,:,t);
    RC{t-1} = [];
    if g.timer, timeElapsed(t-1) = time/npnts; end
end
g.winStartIdx(1) = [];

MODEL.AR = AR;
MODEL.PE = PE;
MODEL.RC = RC;
MODEL.winStartTimes = (g.winStartIdx-1)/EEG.srate;
MODEL.winlen = 1/EEG.srate;
MODEL.winstep = g.downsampleFactor/EEG.srate;
MODEL.morder = g.morder;
MODEL.algorithm = 'kalman';
MODEL.modelclass = 'mvar';
MODEL.timeelapsed = timeElapsed;
MODEL.updatecoeff = g.updatecoeff;
MODEL.updatemode = g.updatemode;
MODEL.downsampleFactor = g.downsampleFactor;
% MODEL.normalize = g.normalize;

