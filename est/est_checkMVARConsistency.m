function stats = est_checkMVARConsistency(varargin)
%
% For a VAR[p] process fit to N time windows, this function returns a
% vector of MVAR model consistency estimates for each time window. The
% percent consistency [1] is an index of the ability for a VAR model, fit
% to data X, to generate data with the same covariance structure as X.
% If Rr and Rs represent the vectorized corss-correlation matrices of the
% real and simulated data, respectively, then the percent consistency is
% given by:
%
% PC = (1 - norm(Rs - Rr)/norm(Rr)) * 100
%
% A value near 100% indicates the model has well-captured the covariance
% structure of the original data
%
% Inputs:
%
%   EEG:            EEG dataset
%   MODEL:          MODEL structure
%   typeproc:       reserved for future use, set to 0 for now
%
% Outputs:
%
%   stats
%       .PC:             Percent Consistency Index [1]
%
% References:
%
% [1] Ding M, Bressler SL, Yang W, Liang H (2000) Short-window spectral
% analysis of cortical event-related potentials by adaptive multivariate
% autoregressive modeling: data preprocessing, model validation, and
% variability assessment. Biol. Cybern. 83:35-45
%
% [2] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
% Theoretical Handbook and User Manual. Chapters 3,6.
% Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% See Also: pop_est_validateMVAR(), est_consistency()
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


g = arg_define([0 2],varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'EEGLAB dataset'), ...
    arg_norep({'MODEL','Model'},mandatory,[],'MVAR MODEL object'), ...
    arg_nogui({'winStartIdx','WindowStartIndices'},[],[],'Starting indices for windows. This is a vector of sample points (start of windows) at which to estimate windowed VAR model. Default is empty (use all windows)','cat','Options'), ...
    arg({'prctWinToSample','WindowSamplePercent'},100,[1 100],'Percent of windows to sample','cat','Options'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );

% commit EEG and MODEL variables to workspace
[data g] = hlp_splitstruct(g,{'EEG','MODEL'});
arg_toworkspace(data);
clear data;


% window size in points
winLenPnts = floor(MODEL.winlen*EEG.srate);

if isempty(g.winStartIdx)
    % starting point of each window (points)
    g.winStartIdx  = round(MODEL.winStartTimes*EEG.srate)+1;
end

if g.prctWinToSample<100
    % randomly select percentage of windows to work with
    randwin = randperm(length(g.winStartIdx));
    randwin = sort(randwin(1:ceil(length(g.winStartIdx)*g.prctWinToSample/100)));
    g.winStartIdx = g.winStartIdx(randwin);
    g.winArrayIndex = randwin;
end

% get the array indices of the windows we are working with
g.winArrayIndex = getindex(MODEL.winStartTimes,(g.winStartIdx-1)/EEG.srate);

if g.verb, h=waitbar(0,sprintf('checking consistency...\nCondition: %s',EEG.condition)); end

numWins = length(g.winStartIdx);
stats.PC = zeros(1,numWins);

for t=1:numWins
    
    % get the array index of the window we are working with
    winArrIdx = g.winArrayIndex(t);
    
    data = squeeze(EEG.CAT.srcdata(:,g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1,:));
    stats.PC(t)= est_consistency(data,MODEL.AR{winArrIdx},MODEL.PE{winArrIdx});
    
    if g.verb,
        waitbar(t/numWins,h,sprintf('checking consistency (%d/%d)...\nCondition: %s',t,numWins,EEG.condition));
    end
end

stats.winStartIdx = g.winStartIdx;
stats.winStartTimes = MODEL.winStartTimes(g.winArrayIndex);
stats.winArrayIndex = g.winArrayIndex;

if g.verb, close(h); end


