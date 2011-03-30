function PC = est_checkMVARConsistency(EEG,MODEL,typeproc,varargin)
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
%   PC:             Percent Consistency Index [1]
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


var = hlp_mergeVarargin(varargin{:});
g = finputcheck(var, hlp_getDefaultArglist('est'), 'est_checkMVARConsistency','ignore');
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

if g.verb, h=waitbar(0,sprintf('checking consistency...\nCondition: %s',EEG.condition)); end

numWins = length(g.winStartIdx);
PC = zeros(1,numWins); 

for t=1:numWins
    
    data = squeeze(EEG.CAT.srcdata(:,g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1,:));
    PC(t)= est_consistency(data,MODEL.AR{t},MODEL.PE{t});

    if g.verb, 
        waitbar(t/numWins,h,sprintf('checking consistency (%d/%d)...\nCondition: %s',t,numWins,EEG.condition));
    end
end

if g.verb, close(h); end


