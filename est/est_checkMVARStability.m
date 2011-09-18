function [stats] = est_checkMVARStability(EEG,MODEL,typeproc,varargin)
%
% Test the stability of a fitted VAR model. See [1-2] for mathematical
% details on testing VAR stability. A stable VAR process is also a
% stationary VAR process [2].
%
% Inputs:
%
%   EEG:        EEGLAB data structure
%   MODEL:      SIFT MODEL structure
%   typeproc:   reserved for future use. Use 0
%
% Optional:
%
%   <Name,Value> pairs containing model fitting parameters. See
%   est_fitMVAR(). Generally, these should be left unspecified.
%
% Outputs:
%
%   stats
%       .stability:  [numwindows x 1] vector of results of stability tests. 1
%                    indicates stable VAR process for that window, 0 indicates
%                    an unstable VAR process.
%
%       .lambda:     [numwindows x nchs*morder] matrix of eigenvalues of VAR
%                    process. All eigenvalues should be < 1 for stable VAR process
%
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual. Chapters 3,6.
%     Available at: http://www.sccn.ucsd.edu/wiki/Sift
% [2] Lutkepohl, H. (2007) New Introduction to Time Series Analysis.
%     Springer.
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD
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

var = hlp_mergeVarargin(varargin{:});
g = finputcheck(var, hlp_getDefaultArglist('est'), 'est_checkMVARStability','ignore','quiet');
if ischar(g), error(g); end


% window size in points
winLenPnts = floor(MODEL.winlen*EEG.srate);

if isempty(g.winStartIdx)
    % starting point of each window (points)
    g.winStartIdx  = floor(MODEL.winStartTimes*EEG.srate)+1;
end

if g.prctWinToSample<100
    % randomly select percentage of windows to work with
    randwin = randperm(length(g.winStartIdx));
    randwin = sort(randwin(1:ceil(length(g.winStartIdx)*g.prctWinToSample/100)));
    g.winStartIdx = g.winStartIdx(randwin);
    MODEL.AR = MODEL.AR(randwin);
    MODEL.PE = MODEL.PE(randwin);
    MODEL.winStartTimes = MODEL.winStartTimes(randwin);
end

if g.verb, h=waitbar(0,sprintf('checking stability...\nCondition: %s',EEG.condition)); end

numWins = length(g.winStartIdx);
stats.stability = zeros(1,numWins);
[nchs Mp] = size(MODEL.AR{1});
stats.lambda = zeros(numWins,Mp);
%lambda = [];
for t=1:numWins
    % rewrite VAR[p] process as VAR[1]
    A = [MODEL.AR{t} ; [eye(nchs*g.morder-nchs,nchs*g.morder-nchs) zeros(nchs*g.morder-nchs,nchs)]];
    stats.lambda(t,:) = log(abs(eig(A)));
    stats.stability(t) = all(stats.lambda(t,:)<0);
    if g.verb,
        waitbar(t/numWins,h,sprintf('checking stability (%d/%d)...\nCondition: %s',t,numWins,EEG.condition));
    end
end

stats.winStartIdx = g.winStartIdx;

if g.verb, close(h); end


