function [IC MODEL params] = est_selModelOrder(EEG,varargin)
%
% Fit a series of MVAR models up to a specified model order and compute the
% model order selection (information) criteria. For additional details see
% [1] and [2].
%
% Inputs:
%
%   EEG                Preprocessed EEG structure. Must contain .CAT
%
% Optional:
%
% 'icselector'         cell array of strings denoting which model order
%                      selection criteria to estimate
%                      'aic': Akaike Information Criterion
%                      'sbc': Swartz Bayes Criterion
%                      'fpe': log of Akaike's Final Prediction Error
%                      'hq': Hannan-Quinn Criterion
% 'algorithm',         string denoting which algorithm to use for model
%                       fitting ('vierra-morf','arfit')
% 'winStartIdx'        vector of sample points (start of windows) at which to estimate windowed VAR model
% 'morder',            [min max] VAR model order to fit
% 'winlen',            window length (sec)
% 'winstep',           window step size (sec)
% 'epochTimeLims',     time range to analyze (sec) where 0 = start of the epoch
% 'prctWinToSample',   percent of time windows to randomly select
% 'verb',              verbosity level (0=no output, 1=text, 2=gui)
% 'normalize'          cell array containing one or more of
%                      {'temporal', 'ensemble'}. This performs ensemble
%                      normalization or temporal normalization (or both)
%                      within each window
%
% Output:
%
%   IC                 a structure containing results of model order selection
%                      IC.selector     - the chosen information criteria
%                      IC.pmin         - the minimum model order tested
%                      IC.pmax         - the maximum model order tested
%                      IC.('sel') contains results for a selector 'sel'.
%                      This consists of subfields
%                           .ic         - [P numwins] matrix of information
%                                         critera for all P model orders tested
%                                         P = morder(2)-morder(1)+1 is the
%                                         number of model orders tested
%                           .minic      - the minimum of ic across model
%                                         orders
%                           .popt       - the model order that minimizes ic
%                           .pelbow     - the model order corresponding to
%                                         the 'elbow' of the information
%                                         criterion (heuristically computed
%                                         as explained in hlp_findElbow())
%                           .elbowic    - the ic value at the 'elbow'
%                           .winStartTimes - the start times of the
%                                            selected windows
%   MODEL               The VAR model fit to pmax.
%   params              The parameters used for model fitting/selection
%
% See Also: pop_est_selModelOrder(), pop_est_fitMVAR(), hlp_findElbow()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapters 3,6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% [2] Lutkepohl, H. (2007) New Introduction to Time Series Analysis.
%   Springer.
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



var = hlp_mergeVarargin(varargin{:});
g = finputcheck(var, [hlp_getDefaultArglist('est'); ...
    {'icselector',          ''     {}   {'sbc','aic','fpe','hq','ris'};}], ...
    'est_selModelOrder','ignore','quiet');
if ischar(g), error(g); end
if nargout > 2, params = g; end

if length(g.morder)<2
    error('''morder'' should contain a minimum and maximum model order');
end

pmin        = g.morder(1);
pmax        = g.morder(2);
npnts       = EEG.trials*floor(g.winlen*EEG.srate);
nbchan      = EEG.CAT.nbchan;

if ~isempty(g.icselector) && ischar(g.icselector)
    g.icselector = {g.icselector};
end

if g.verb
    est_dispMVARParamCheck(EEG,g);
end

if ismember(g.algorithm,{'vieira-morf-cpp','arfit'})
    % for these methods, we have to fit a separate MODEL for each model order
    for p=pmin:pmax
        VARtmp(p-pmin+1) = est_fitMVAR(EEG,0,g,'morder',p);
    end
    numWins         = length(VARtmp(1).winStartTimes);
    winStartTimes   = VARtmp(1).winStartTimes;
    for t=1:numWins
        for p=pmin:pmax,
            % extract noise covariance matrix for each model order and window
            MODEL.PE{t}(:,p*nbchan+(1:nbchan)) = VARtmp(p-pmin+1).PE{t}(end-nbchan+1:end,1:nbchan);
        end
    end;
else
    % fit MVAR model up to maximum model order
    MODEL           = est_fitMVAR(EEG,0,g,'morder',pmax);
    numWins         = length(MODEL.winStartTimes);
    winStartTimes   = MODEL.winStartTimes;
end

% initialize some variables
[sbc fpe aic hq ris]    = deal(nan*ones(pmax-pmin+1,numWins));
nparams = nbchan^2.*(pmin:pmax);

for t=1:numWins
    
    % CALCULATE INFORMATION CRITERIA
    
    ne = npnts-(pmin:pmax);
    logdp = zeros(1,pmax-pmin+1);
    
    for p=pmin:pmax,
        % Get logarithm of determinant for each model order
        logdp(p-pmin+1) = log(det(MODEL.PE{t}(:,p*nbchan+(1:nbchan))*(npnts-p)));
    end;
    
    
    % Schwarz's Bayesian Criterion
    sbc(:,t) = logdp + (log(ne).*nparams./ne);   % TM
    
    % Akaike Information Criterion
    aic(:,t) = logdp + 2.*nparams./ne;   % TM
    
    % logarithm of Akaike's Final Prediction Error
    fpe(:,t) = logdp + nbchan*log(ne.*(ne+nbchan*(pmin:pmax))./(ne-nbchan*(pmin:pmax)));  % TM
    
    % Hannan-Quinn criterion
    hq(:,t) = logdp + nparams.*2.*log(log(ne))./ne;  % TM
    
    % Rissanen criterion
    ris(:,t) = logdp + (nparams./ne).*log(ne);
    
    for i=1:length(g.icselector)
        % get index iopt of order that minimizes the order selection
        % criterion specified in g.icselector
        sel = g.icselector{i};
        ic = eval([sel '(:,t);']);
        [minic.(sel)(t) iopt] = min(ic);
        popt.(sel)(t) = pmin + iopt-1; % estimated optimum order
        
        
        % get model order corresponding to the "elbow" of the order
        % selection criterion. An "elbow" is found using a geometric
        % heuristic (see hlp_findElbow() for details)
        [elbowic.(sel)(t) iopt] = hlp_findElbow(ic);
        pelbow.(sel)(t) = pmin + iopt-1; % estimated optimum order
    end
    
end

% store the information criteria in output structure
for i=1:length(g.icselector)
    sel = g.icselector{i};
    eval(['IC.(sel).ic = ' sel ';']);
    IC.(sel).minic = minic.(sel);
    IC.(sel).popt = popt.(sel);
    IC.(sel).pelbow = pelbow.(sel);
    IC.(sel).elbowic = elbowic.(sel);
end

IC.selector = g.icselector;
IC.pmin = pmin;
IC.pmax = pmax;
IC.winStartTimes = winStartTimes;


