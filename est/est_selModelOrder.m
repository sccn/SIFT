function [IC g] = est_selModelOrder(EEG,varargin)
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
%   cfg              The parameters used for model fitting/selection
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


g = arg_define([0 1],varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'EEGLAB dataset'), ...
    arg_subswitch({'modelingApproach','ModelingApproach'},'Segmentation VAR', ...
    { ...
    {'Segmentation VAR' @est_fitMVAR}, ...   % {'Kalman Filtering' {arg('arg_dummy',[],[],'replace this with @est_fitMVARKalman')}}, ...
    }, ...
    'Select a modeling approach and define parameters. Make sure to use the same parameters you intend to use for the final modeling step.', ...
    'cat','Modeling Parameters', ...
    'suppress',{'verb','timer','ModelOrder','WindowSamplePercent','EpochTimeLimits','WindowStartIndices'}), ...
    arg({'morderRange','ModelOrderRange','modelOrderRange'},[1 30],[],'VAR model order range.','cat','Modeling Parameters','shape','row'), ...
    arg({'downdate','Downdate'},true,[],'Downdate model. If selected (and modeling approach supports it), a model of the highest desired order will be fit and then lower-order models will be approximated by downdating. This can be much faster than successively fitting models up to the maximum order. Note however, (1) this only renders an approximation and (2) this requires that there is sufficient data to adequately fit a model of the highest order.', ...
    'cat','Modeling Parameters'), ...
    arg({'icselector','InformationCriteria'},{'sbc','aic','fpe','hq','ris'},{'sbc','aic','fpe','hq','ris'},'Order selection criteria. This specifies the information criteria to use for order selection.','cat','Modeling Parameters','type','logical'), ...
    arg_nogui({'winStartIdx','WindowStartIndices'},[],[],'Starting indices for windows. This is a vector of sample points (start of windows) at which to estimate windowed VAR model','cat','Data Selection'), ...
    arg_nogui({'epochTimeLims','EpochTimeLimits'},[],[],'Epoch time limits (sec). This is relative to event time (e.g. [-1 2]). Default is the full epoch time range','cat','Data Selection'), ...
    arg({'prctWinToSample','WindowSamplePercent'},100,[1 100],'Percent of windows to sample','cat','Data Selection'), ...
    arg_subtoggle({'plot','PlotResults'},{},@vis_plotOrderCriteria,'Plot results','suppress',{'InformationCriteriaToPlot','TitleString'}), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );

% commit EEG variable to workspace
[data g] = hlp_splitstruct(g,{'EEG'});
arg_toworkspace(data);
clear data;

% do some error checking
if length(g.morderRange)<2
    error('''morderRange'' should contain a minimum and maximum model order');
end
if g.morderRange(1)<1
    error('minimum allowable model order is 1');
end

pmin        = g.morderRange(1);
pmax        = g.morderRange(2);
nbchan      = EEG.CAT.nbchan;

if ~isempty(g.icselector) && ischar(g.icselector)
    g.icselector = {g.icselector};
end

% if g.verb
%     est_dispMVARParamCheck(EEG,g.modelingApproach,true);
% end

% determine whether we can downdate the model
if g.downdate
    if strcmp(g.modelingApproach.arg_selection,'Kalman Filtering') ...
            || (strcmp(g.modelingApproach.arg_selection,'Segmentation VAR') ...
            && ~any(strcmpi(g.modelingApproach.algorithm.arg_selection,{'vieira-morf'})))
        % the modeling approach does not support downdating
        fprintf(['WARNING: The selected modeling approach/algorithm does not support downdating\n' ...
            'Switching to exhaustive search mode\n']);
        g.downdate = false;
    end
end

% % determine which windows to use
% if isempty(g.winStartIdx)
%     % starting point of each window (points)
%     g.winStartIdx  = round(MODEL.winStartTimes*EEG.srate)+1;
% end
% 
% % randomly select percentage of windows to work with
% if g.prctWinToSample<100
%     randwin = randperm(length(g.winStartIdx));
%     randwin = sort(randwin(1:ceil(length(g.winStartIdx)*g.prctWinToSample/100)));
%     g.winStartIdx = g.winStartIdx(randwin);
%     g.prctWinToSample = 100;
% end


% get the handle of the modeling function
switch g.modelingApproach.arg_selection
    case 'Segmentation VAR'
        modelingFuncName = 'est_fitMVAR';
    case 'Kalman Filtering'
        modelingFuncName = 'est_fitMVARKalman';
end
    
if ~g.downdate
    % sequentially fit a separate MODEL for each model order
    
    if g.verb==2
        h = waitbar(0,'Sequentially searching model order range...');
    end
    
    cnt = 0; tot = (pmax-pmin)+1;
    for p=pmin:pmax
        
        % initialize random number generator with same seed
        % so that we get the same random sequence of windows for 
        % each model order
        if g.prctWinToSample<100
            rng(1);  
        end
        
        cnt = cnt + 1;
        
        VARtmp(p-pmin+1) = feval(modelingFuncName,EEG,g.modelingApproach, ...
                                'ModelOrder',p,'verb',g.verb, ...
                                'epochTimeLims',g.epochTimeLims, ...
                                'winStartIdx',g.winStartIdx, ...
                                'prctWinToSample',g.prctWinToSample);
        
        if g.verb==2
            waitbar(cnt/tot,h, ...
                sprintf('Sequentially searching model order range (%d/%d)...',cnt,tot));
        end
    end
    
    % cleanup
    if g.verb==2
        delete(h); 
    end
    
    numWins         = length(VARtmp(1).winStartTimes);
    winStartTimes   = VARtmp(1).winStartTimes;
    winlen          = VARtmp(1).winlen;
    
    % extract noise covariance matrix for each model order and window
    for t=1:numWins
        for p=pmin:pmax,
            MODEL.PE{t}(:,p*nbchan+(1:nbchan)) = VARtmp(p-pmin+1).PE{t}(end-nbchan+1:end,1:nbchan);
        end
    end;
    MODEL.winlen = VARtmp(end).winlen;
else
    % fit MVAR model up to maximum model order 
    % Lower orders will be approximated by downdating
    MODEL           = feval(modelingFuncName,EEG,g.modelingApproach, ...
                                'ModelOrder',pmax,'verb',g.verb, ...
                                'epochTimeLims',g.epochTimeLims, ...
                                'winStartIdx',g.winStartIdx, ...
                                'prctWinToSample',g.prctWinToSample);
    numWins         = length(MODEL.winStartTimes);
    winStartTimes   = MODEL.winStartTimes;
    winlen          = MODEL.winlen;
end

% initialize some variables
[sbc fpe aic hq ris]    = deal(nan*ones(pmax-pmin+1,numWins));
nparams = nbchan^2.*(pmin:pmax);

npnts       = EEG.trials*max(1,floor(winlen*EEG.srate));

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

IC.modelFitting.modelingFuncName = modelingFuncName;
IC.modelFitting.modelingArguments = g.modelingApproach;
IC.selector = g.icselector;
IC.pmin = pmin;
IC.pmax = pmax;
IC.winStartTimes = winStartTimes;


% plot results (if desired)
if g.plot.arg_selection
    g.plot.conditions = EEG.condition;
    vis_plotOrderCriteria(IC,g.plot);
end
