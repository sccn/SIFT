function [AR PE out] = mvar_scsa_em(varargin)

% Algorithm: Group Lasso (ADMM)
%
% Description:
%
% This method infers a sparsely connected multi-
% variate autoregressive (VAR) model under a 
% gaussian noise assumption.
%
% VAR[p] coefficients inferred using Group Lasso 
% (L1,2 regularization) via the Alternating 
% Direction Method of Multipliers [1]. 
% The p VAR coefficients describing interactions 
% between two processes at time lags [t-1...t-p]
% are grouped together using an L2 norm which 
% penalizes large coefficients. An L1 penalty is
% then applied to the coefficient groups. This
% jointly prunes entire groups of coefficients
% by setting the entire group to zero. The result
% is a connectivity graph with sparse structure
% (most processes are strictly non-interacting)
% while regularizing (smoothing) the coefficient
% sequences describing surviving non-zero 
% interactions.
% These constraints allow us to uniquely solve 
% highly under-determined systems (e.g. many 
% more parameters than observations).
%
% If Y = D(y,p) is a p-lag delay embedding of multi-
% variate data vector y(t), X is a sparse block-
% diagonal matrix of lagged copies of the delay 
% embedded data, and A is an augmented matrix of
% VAR coefficients, then we may adopt the structural 
% model:
%
% Y = XA + e,  for gaussian noise e
%
% We then seek to solve the optimization problem
%
% A_hat = argmin_A{0.5||Y-XA||_2^2 + L*sum(||A_i||_2)}
%
% where A_i contains the i^th group of AR parameters 
% e.g. the set of VAR weights {A(i,j)} describing auto-
% regression (conditional linear dependence) of X_i 
% onto X_j. L is the regularization parameter.
%
% This algorithm is a good choice if you have few 
% data samples and many channels/sources and/or a 
% high model order.
%
% Author Credits:
%
% This implementation is based on Matlab ADMM
% examples from Stephen Boyd's website [2]
%
% References and Code:
%
% [1] Boyd, Parikh, Chu, Pelato and Eckstein et al. Foundations and Trends in Machine Learning 3(1):1-122,2011
% [2] http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
% Dependencies: admm_gl()
%
% ------------------------------------------------------------------------
% INPUTS:
%   data:       the data (nchs x npnts)
%   p:          the model order
%   lambda:     regularization parameter
%   rho:        the augmented Lagrangian parameter
%   alpha:      the over-relaxation parameter (typical 
%               values for alpha are between 1.0 and 1.8)
%   AR0:        initial solution (default: zeros)
%   verb:       true/false - verbosity
% OUTPUTS:
%   AR:         [nchs x nchs*p] VAR coefficient matrix
%   PE:         [nchs x nchs] estimated noise covariance matrix
%
% See Also: est_fitMVAR(), mvar_dalSCSA(), mvar_ridge(), mvar_vieiramorf(),
%           est_fitMVARKalman()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift/
% [2] S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein, "Distributed 
%       Optimization and Statistical Learning via the Alternating Direction 
%       Method of Multipliers". Foundations and Trends in Machine Learning, 
%       Michael Jordan, Editor in Chief, 3(1):1-122, 2011.
%       http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
% [3] http://www.stanford.edu/~boyd/papers/admm/group_lasso/group_lasso_example.html
%
% Author: Tim Mullen 2012, SCCN/INC, UCSD. 
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

persistent initAR;

% g = arg_define([0 1],varargin, ...
%                 arg_norep({'data','Data'},mandatory,[],'Data Matrix. Dimensions are [nchs x npnts].'), ...
%                 arg({'morder','ModelOrder','p'},10,[],'VAR Model order'), ...
%                 arg_subtoggle({'warmStart','WarmStart'},[], ...
%                 {...
%                     arg({'initState','InitialState'},[],[],'Initial ADMM state. This is a structure with fields ''z_init'' and ''u_init'', which represent, respectively, the initial state vector and dual vector and are both of dimension [morder*(nchs^2) x 1]. If empty, the first state is initialized to zeros.') ...
%                 },'Warm start. The previously estimated state will be used as a starting estimate for successive operations. Alternately, you may provide an non-empty initial state structure via the ''InitialState'' argument.'), ...
%                 arg({'normcols','NormCols'},'none',{'none','norm','zscore'},'Normalize columns of dictionary'), ...
%                 arg({'groupDiags','GroupAutoConnections','GroupDiags'},false,[],'Group auto-connections. All auto-connections for all channels will be penalized jointly as a separate group. Note, this can slow down model estimation as the design matrix will not longer be diagonal block-toeplitz so we cannot (yet) exploit block-redundancy in the design matrix'), ...
%                 arg_norep({'AR0','InitialState'},[],[],'DEPRECATED. Initial VAR coefficient matrix','shape','matrix','type','expression'), ...
%                 arg_sub({'admm_args','ADMM_Options'},[],@admm_gl,'Options for ADMM algorithm') ...
%                 );
%                 
% arg_toworkspace(g);

g = varargin{1};
  
[nchs npnts ntr] = size(g.data);
p = g.morder;

% regularization constant
lambda = g.lambda;

% initialize state
if ~isempty(g.warmStart)
    initAR.B = g.warmStart.B;
    initAR.H = g.warmStart.H;
    runCSA = 0;
else
    % reset initAR
    initAR.B = eye(nchs);
    initAR.H = zeros(nchs, nchs, p);
    runCSA = 1;
end

% perform dimensionality reduction to g.PCA largest principal components
if ~isempty(g.PCA) 
  % center each channel/trial
  g.data = bsxfun(@minus, g.data, sum(g.data,2)/npnts);
  
  % indicate that the data have been made zero-mean before running PCA and
  % SCSA. The same must be done if one wants to use out.pcafilt or 
  % out.scsafilt to recover the sources
  out.zeromean = 1;
  
  if isstruct(g.PCA)
    out.pcafilt = g.PCA.pcafilt;
    out.pcapat  = g.PCA.pcapat;
    g.PCA       = size(out.pcafilt, 1);
    g.data      = reshape(out.pcafilt*g.data(:,:), [g.PCA, npnts, ntr]);
  elseif g.PCA > 0 && g.PCA < nchs
    % perform PCA
    % first compute average data covariance matrix across trials...
    % ...concat trials (reshape to [nchs x npnts*ntr])
    cdata = g.data(:,:);
    % ...compute (unbiased) mean covariance
    C     = cdata*cdata'/(ntr*(npnts-1));
    % compute eigendecomposition
    [V D] = eig(C);
    d     = diag(D);
    out.varianceexplained = sum(d(end-g.PCA+1:end))/sum(d);
    out.pcafilt           = diag(1./sqrt(d(end-g.PCA+1:end)))*V(:, end-g.PCA+1:end)';
    out.pcapat            = V(:, end-g.PCA+1:end)*diag(sqrt(d(end-g.PCA+1:end)));
    % store principal components
    g.data                = reshape(out.pcafilt*cdata, [g.PCA, npnts, ntr]);
    clear cdata;
  end
end

opts = struct('start', initAR, 'lambda', lambda, 'runCSA', runCSA, 'emstop', 1e-6);
[out.scsafilt, AR, out.f, out.S] = scsa_em(g.data, p, opts);

out.scsapat = inv(out.scsafilt);

if ~isempty(g.PCA) && g.PCA > 0 && g.PCA < nchs
  out.pat = out.pcapat/out.scsafilt;
  out.filt = out.scsafilt*out.pcafilt;
else
  out.pat = out.scsapat;
  out.filt = out.scsafilt;
end

%% estimate noise covariance matrix
if nargout>1
    res = est_mvarResiduals(out.S,AR,zeros(1,g.PCA));
    res = res(:,:);
    PE = cov(res',1);
end


