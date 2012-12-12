function [AR PE] = mvar_glADMM(varargin)

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
%       Michael Jordan, Editor in Chief, 3(1):1?122, 2011.
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

g = arg_define([0 1],varargin, ...
                arg_norep({'data','Data'},mandatory,[],'Data Matrix. Dimensions are [nchs x npnts].'), ...
                arg({'morder','ModelOrder'},10,[],'VAR Model order'), ...
                arg_subtoggle({'warmStart','WarmStart'},[], ...
                {...
                    arg({'initState','InitialState'},[],[],'Initial ADMM state vector. The dimension is [morder*(nchs^2) x 1]. If empty, the first state is initialized to zeros.') ...
                },'Warm start. The previously estimated state will be used as a starting estimate for successive operations. Alternately, you may provide an non-empty initial state structure via the ''InitialState'' argument.'), ...
                arg({'normcols','NormCols'},'none',{'none','norm','zscore'},'Normalize columns of dictionary'), ...
                arg({'groupDiags','GroupAutoConnections','GroupDiags'},false,[],'Group auto-connections. All auto-connections for all channels will be penalized jointly as a separate group. Note, this can slow down model estimation as the design matrix will not longer be diagonal block-toeplitz so we cannot (yet) exploit block-redundancy in the design matrix'), ...
                arg_norep({'AR0','InitialState'},[],[],'DEPRECATED. Initial VAR coefficient matrix','shape','matrix','type','expression'), ...
                arg_sub({'admm_args','ADMM_Options'},[],@admm_gl,'Options for ADMM algorithm') ...
                );
                
% arg_toworkspace(g);
  
[nchs npnts ntr] = size(g.data);
p = g.morder;

% initialize state
if g.warmStart.arg_selection && ~isempty(g.warmStart.initState)
    initAR = g.warmStart.initState;
elseif ~g.warmStart.arg_selection
    % reset initAR
    initAR = zeros(p*nchs^2,1);
end
if size(initAR,1) ~= p*nchs^2
    % dimensions have changed, reset state
    fprintf('mvar_glADMM: model dimensions changed -- resetting state\n');
    initAR = zeros(p*nchs^2,1);
end

% Indices for the diagonal elements (self connection)
diagIdx = vec(repmat((1:p*(nchs+1):p*nchs^2), p, 1) + repmat((0:(p-1))', 1, nchs))';
% Indices for the off-diagonal elements (connection to others)
offDiagIdx = setdiff(1:p*nchs^2, diagIdx);

% Reshape the design matrix into VAR[1] format
%
% X will be of size [ntr*(npnts-p) x nchs*p]
% Y will be of size [ntr*(npnts-p) x nchs]
%
% X consists of ntr vertically stacked Toeplitz matrices.
% Suppose ntr==1, M=nchs. Then X consists of M blocks [X1,...,XM]
% where Xi is a matrix of size [npnts-p x p] containing 
% the (1...p) order delay embedding of the data for channel i.
% That is, Xi(:,k) is the data vector for channel i delayed by p+1-k 
% samples (smaller k => more delay)
%
% Note that Y = X0.
blocklen   = npnts-p;
blockwidth = p*nchs;
X = zeros(blocklen,blockwidth,ntr);
Y = zeros(blocklen,nchs,ntr);
for itr = 1:ntr
    % initialize delay-embedding blocks
    Xi = zeros(blocklen,nchs,p);
    for jj = 1:p
        Xi(:,:,jj) = squeeze(g.data(:, p+1-jj:end-jj, itr))';
    end
    % permute to [npnts x p x nchs] so each page is a 
    % delay-embedding block for a given channel...
    Xi = permute(Xi, [1 3 2]); 
    % ... and concatenate blocks horizontally to form 2D design mat
    % for this trial X(:,:,itr) = [X1 X2 ... XM]
    X(:,:,itr) = reshape(Xi, blocklen, blockwidth);
    
    % extract 0-lag data matrix for this trial...
    Y(:,:,itr) = g.data(:, p+1:end,itr)';
end
clear Xi;
% reshape X and Y to stack trials vertically and form final 2D matrix
X = permute(X,[1 3 2]);  % [npnts x trials x nchs]
Y = permute(Y,[1 3 2]);
X = reshape(X,blocklen*ntr,blockwidth);
Y = reshape(Y,blocklen*ntr,nchs);

% %% OLD METHOD
% % Reshape the design matrix into VAR[1]
% X = [];
% Y = [];
% for itr = 1:ntr
%     Xi = [];
%     for jj = 1:p
%         Xi = cat(3, Xi, squeeze(g.data(:, p+1-jj:end-jj, itr))');
%     end
%     X = cat(1, X, reshape(permute(Xi, [1 3 2]), npnts-p, nchs*p));
%     Y = cat(1,Y,g.data(:, p+1:end,itr)');
% end

% target vector
Y = vec(Y);

% normalization
switch lower(g.normcols)
    case 'zscore'
        % standardize columns of X (unit variance)
        X = zscore(X,1,1);
        Y = zscore(Y);
    case 'norm'
        % normalize columns of X (unit norm)
        X     = bsxfun(@rdivide,X,sqrt(sum(X.^2)));
        Y     = Y./norm(Y);
end

% predictor design matrix
X = blkdiageye(sparse(X),nchs);
if g.groupDiags
    X = [X(:, offDiagIdx) X(:, diagIdx)];
    blks = [p*ones(1,nchs*(nchs-1)),p*nchs];
else
    blks = p*ones(1,nchs^2);
end


% Apply the ADMM method for group lasso estimation:
% group penalize AR coefficients with different time-lags
% between sources (nchs*(nchs-1) groups of size p). 
[initAR] = admm_gl('A',X,  ...
                   'y',Y, ...
                   'blks',blks, ...
                   g.admm_args, ...
                   'z_init',initAR, ...
                   'designMatrixBlockSize',fastif(g.groupDiags,[],[blocklen*ntr blockwidth]));

% assemble coefficient matrices
AR = zeros(p,nchs,nchs);
if g.groupDiags
    AR(offDiagIdx) = vec(initAR(1:(p*nchs*(nchs-1))));          % non-diagonal elements
    AR(diagIdx) = vec(initAR(((p*nchs*(nchs-1))+1):end));       % diagonal elements
else
    AR(:) = initAR;
end

AR = permute(full(AR), [3 2 1]);   % the recovered (nchs x nchs x p) connectivity matrix
AR = reshape(AR,[nchs nchs*p]);

%% estimate noise covariance matrix
if nargout>1
    res = est_mvarResiduals(g.data,AR,zeros(1,nchs));
    res = res(:,:);
    PE = cov(res');
end

% if nargout>2
%     argsout.initAR = initAR;
% end

function v = vec(x)
v = x(:);

