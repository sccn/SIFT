
function [stats acv params] = est_checkMVARWhiteness(EEG,MODEL,typeproc,varargin)

% Tests for whiteness of residuals of VAR model. Residuals are 'white' if
% they are statistically uncorrelated (e.g. a white noise process). White
% residuals indicate the second-order dynamics of the data is adequately
% captured by the VAR model.
%
% If stats.pval > 0.05, we cannot reject the null hypothesis (that
% the data is white) at the 0.05 significance level, in other
% words the probability of mistakenly assuming the data is white is less
% than 5%, or, equivalently, there is at least a 95% probability that we 
% are correct in assuming the data is white.
%
% INPUT:
%
%   EEG:         EEG data structure
%   MODEL:       MODEL structure returned by est_fitMVAR or similar
%
% OPTIONAL:
%   'prct2sample':      {def: 100} percent of windows to sample [0 100]
%   'alpha':            {def: 0.05} whiteness significance level
%   'whitenessCriteria': Whiteness tests. Can be one or more of 
%                       'ACF','Ljung-Box','Box-Pierce','Li-McLeod'
%   'verb':             {def: true} enable verbosity
%   'numLagsAcf':       {def: 50} number of autocorrelation lags for 
%                       whiteness tests
% OUTPUT:
%   stats:
%    .[whitenessCriteria{i}]
%       .w:         vector of logicals. w(t)=1 if residuals for window t are 
%                   white, 0 otherwise
%       .pval:      vector of whitness p-values. If pval(t) < alpha,
%                   residuals are white at the alpha-level for window t.
%                   Alpha should be in [0 1]
%       .value:     values of the statistic
%       .fullname:  original name of test statistic
%       .winStartTime: start times of each window (sec)
%   acv:     [nchs x nchs x 2*numLagsAcf+1 x numwins] matrix of autocorrelation
%            coeffs for lags 0:numLagsAcf
%   params:  struct of options used
%
%
% See Also: est_fitMVAR(), pop_est_validateMVAR()
%
% References:
%
% [1] Li, W. K., and A. I. McLeod, (1981). Distribution of the
%     Residual Autocorrelations in Multivariate ARMA Time Series
%     Models, J. Roy. Stat. Soc. B, 43, 231--239.      
% [2] Lutkepohl, H. (2007) New Introduction to Time Series Analysis.
%     Springer.
% [3] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual. Chapters 3,6. 
%     Available at: http://www.sccn.ucsd.edu/wiki/Sift
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
myargs = {'whitenessCriteria'   'cell'  {'ljungbox','acf','boxpierce','lmp'}   {'lmp'}};
g = finputcheck(var, [myargs ; hlp_getDefaultArglist('est')], 'est_checkWhiteness','ignore','quiet');
if ischar(g), error(g); end
if nargout > 2, params = g; end
g.whitenessCriteria = lower(g.whitenessCriteria);
if isempty(g.morder)
    g.morder = MODEL.morder;
end

% initialize vars
acv = [];
if ~isfield(g,'alpha') || isempty(g.alpha)
    g.alpha = 0.05;  % significance threshold
end

% window size in points
winLenPnts = floor(g.winlen*EEG.srate); 
[nchs npnts ntr] = size(EEG.CAT.srcdata);
% total number of points used to estimate model
npntsTotal = ntr*(winLenPnts-g.morder);            

if isempty(g.winStartIdx)
    % starting point of each window (points)
    g.winStartIdx  = round(MODEL.winStartTimes*EEG.srate)+1;    
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
numWins = length(g.winStartIdx);
statcorrection = 1; %nchs;   % bonferroni correction for multiple comparisons

% number of lags for acf whiteness testing
if ~isfield(g,'numLagsAcf') || isempty(g.numLagsAcf)
    numLagsAcf = 50; min(winLenPnts-g.morder-1, max([50 6*g.morder])); %max(30,g.morder+1); %g.morder-1; %min(max(20,g.morder-1), winLenPnts-1); %min(5,g.morder);  % g.morder? 
end

if g.verb, fprintf('Using %d lags for covariance estimation.\n',numLagsAcf); end
if g.verb, h=waitbar(0,sprintf('checking whiteness...\nCondition: %s',EEG.condition)); end


for t=1:numWins
    
    if any(ismember(g.whitenessCriteria,{'acf','ljungbox','boxpierce','limcleod'}))
        % calculate residuals
        residuals = est_mvarResiduals(squeeze(EEG.CAT.srcdata(:,g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1,:)), ...
                                      EEG.CAT.MODEL.AR{t},zeros(1,EEG.CAT.nbchan));
    end
    
    if any(ismember(g.whitenessCriteria,{'ljungbox','acf','boxpierce','limcleod'}))

%         if size(MODEL.PE{t},2)>nchs
%             C = MODEL.PE{t}(:,nchs*g.morder+1:nchs*(g.morder+1));
%         else
%             C = MODEL.PE{t};
%         end
% 
%         if g.morder>numLagsAcf
%             
%             % THIS GETS THE COVARIANCE MATRIX OF THE VAR PROCESS -- NOT THE
%             % RESIDUALS...
%             
%             acvf2{t} = est_calcInvCovMat(MODEL.AR{t},C,false);   % get autocovariance matrix
%             acvf2{t} = acvf2{t}(1:nchs,:);    % extract positive lags (acvf is block-toeplitz and hermitian so only need top block-row)
%             % calculate autocorrelation matrix
%             acf2{t} = zeros(nchs,nchs*(numLagsAcf+1));
%             D = inverse(sqrt(diag(diag(acvf2{t}(1:nchs,1:nchs)))));  %\eye(nchs);
%             for k=0:numLagsAcf
%                 acf2{t}(:,k*nchs+1:(k+1)*nchs) = D*acvf2{t}(1:nchs,(1:nchs)+k*nchs)*D;
%             end
%             
%         end
            
        % calculate correlation matrix for each trial
        [acvf{t} acf{t}] = deal(zeros(ntr,nchs,nchs*(numLagsAcf+1)));
        for tr=1:ntr
            % R = [Rs1s1 Rs1s2 Rs1s3 Rs2s1 Rs2s2 Rs2s3 Rs3s1 Rs3s2 Rs3s3]
            tmpacvf = xcov(squeeze(residuals(:,:,tr))',numLagsAcf,'biased');
            tmpacvf = tmpacvf(numLagsAcf+1:end,:);         % extract positive lags

            for lag = 1:numLagsAcf+1
                acvf{t}(tr,:,(1:nchs)+nchs*(lag-1)) = reshape(squeeze(tmpacvf(lag,:)),[nchs nchs])';
            end

            % compute autocorrelation function
            D = diag(1./sqrt(diag(squeeze(acvf{t}(tr,1:nchs,1:nchs)))));
            for lag = 1:numLagsAcf+1
                acf{t}(tr,:,(1:nchs)+nchs*(lag-1)) = D*squeeze(acvf{t}(tr,:,(1:nchs)+nchs*(lag-1)))*D;
            end

        end

        acf{t} = squeeze(mean(acf{t},1));
        acvf{t} = squeeze(mean(acvf{t},1));
            
    end
    
    % compute whiteness criteria
    % TODO: performance can be improved by avoiding recomputation of
    % similar statistics 
    for i=1:length(g.whitenessCriteria)
        switch lower(g.whitenessCriteria{i})
            
            case 'acf'
                % test residual autocorrelations
                sigthresh = 1.96/sqrt(npntsTotal);    % 0.95 confidence interval for acvf
                numNonSigCoeffs = sum(sum(abs(acf{t}(:,nchs+1:end))<sigthresh));    % count how many coefficients are within the bounds for 'white' residuals
                numCoeffs = numel(acf{t}(:,nchs+1:end));                            % count the total number of coefficients
                stats.acf.pval(t) = numNonSigCoeffs / numCoeffs;                    % estimate the probability for a coefficient to be inside the bounds (probability of whiteness)
                stats.acf.w(t) = stats.acf.pval(t) > 1-g.alpha;                       
                stats.acf.acffun{t} = acf{t};
                stats.acf.fullname = 'ACF';
                stats.acf.winStartTimes = g.winStartIdx*EEG.srate;
            case 'ljungbox' 
                % ljung-box (modified portmanteau) test for residual autocorrelation up to lag h
                Qh=0;
                C0inv = inverse(acvf{t}(1:nchs,1:nchs));
                for k=1:numLagsAcf
                    Qh = Qh + 1/(npntsTotal-k) * trace(acvf{t}(1:nchs,(1:nchs)+k*nchs)'*C0inv*acvf{t}(1:nchs,(1:nchs)+k*nchs)*C0inv);
                end
                stats.ljungbox.value(t) = Qh*npntsTotal*(npntsTotal+2);
                stats.ljungbox.pval(t) = 1-chi2cdf(stats.ljungbox.value(t),(nchs^2)*(numLagsAcf-g.morder));
                stats.ljungbox.w(t)= stats.ljungbox.pval(t)>g.alpha;
                stats.ljungbox.fullname = 'Ljung-Box';
                stats.ljungbox.winStartTimes = g.winStartIdx*EEG.srate;
            case 'boxpierce' 
                % box-pierce portmanteau test for residual autocorrelation up to lag h
                Qh=0;
                C0inv = inverse(acvf{t}(1:nchs,1:nchs));
                for k=1:numLagsAcf
                    Qh = Qh + trace(acvf{t}(1:nchs,(1:nchs)+k*nchs)'*C0inv*acvf{t}(1:nchs,(1:nchs)+k*nchs)*C0inv);
                end
                stats.boxpierce.value(t) = Qh*npntsTotal;   % Qh*npntsTotal^2  is the modified portmanteau test (Lutkepohl, p.171)
                stats.boxpierce.pval(t) = 1-chi2cdf(stats.boxpierce.value(t),(nchs^2)*(numLagsAcf-g.morder));
                stats.boxpierce.w(t)= stats.boxpierce.pval(t)>g.alpha;
                stats.boxpierce.fullname = 'Box-Pierce';
                stats.boxpierce.winStartTimes = g.winStartIdx*EEG.srate;
            case 'limcleod'
                % Li-McLeod portmanteau test for residual autocorrelation up to lag h
                Qh=0;
                C0inv = inverse(acvf{t}(1:nchs,1:nchs));
                for k=1:numLagsAcf
                    Qh = Qh + trace(acvf{t}(1:nchs,(1:nchs)+k*nchs)'*C0inv*acvf{t}(1:nchs,(1:nchs)+k*nchs)*C0inv);
                end
                stats.limcleod.value(t) = Qh*npntsTotal + nchs^2*numLagsAcf*(numLagsAcf+1)/(2*npntsTotal);
                stats.limcleod.pval(:,t) = 1-chi2cdf(stats.limcleod.value(t),(nchs^2)*(numLagsAcf-g.morder));
                stats.limcleod.w(t) = stats.limcleod.pval(t) > g.alpha;    % using bonferroni correction for multiple tests (across channels)
                stats.limcleod.fullname = 'Li-McLeod';
                stats.limcleod.winStartTimes = g.winStartIdx*EEG.srate;
            otherwise
                fprintf('Unknown whiteness test %s',g.whitenessCriteria{i}); 
        end
    end
    
    if g.verb
        waitbar(t/numWins,h,sprintf('checking whiteness (%d/%d)...\nCondition: %s',t,numWins,EEG.condition)); 
    end
                
end

stats.winStartIdx = g.winStartIdx;
stats.alpha = g.alpha;

if g.verb, close(h); end

