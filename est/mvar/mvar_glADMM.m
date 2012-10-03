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

verb = arg_extract(varargin,{'verb','Verbosity'},[],0);

g = arg_define([0 1],varargin, ...
                arg_norep({'data','Data'},mandatory,[],'Data Matrix. Dimensions are [nchs x npnts].'), ...
                arg({'morder','ModelOrder'},10,[],'VAR Model order'), ...
                arg_nogui({'AR0','InitialState'},[],[],'Initial VAR coefficient matrix','shape','matrix','type','expression'), ...
                arg({'verb','Verbosity'},verb,[],'Verbose output'), ...
                arg_sub({'admm_args','ADMM_Options'},[],@admm_gl,'Options for ADMM algorithm') ...
                );
                
% arg_toworkspace(g);
  
[nchs npnts ntr] = size(g.data);
p = g.morder;

% initialize state
if ~isempty(g.AR0)
    initAR = g.AR0;
elseif isempty(initAR)
    initAR = zeros(p*nchs^2,1);
end
if size(initAR,1) ~= p*nchs^2
    % dimensions have changed, reset state
%     if g.verb
        fprintf('mvar_glADMM: model dimensions changed -- resetting state\n');
%     end
    initAR = zeros(p*nchs^2,1);
end

% Reshape the design matrix into VAR[1]
X = [];
Y = [];
for itr = 1:ntr
    Xi = [];
    for jj = 1:p
        Xi = cat(3, Xi, squeeze(g.data(:, p+1-jj:end-jj, itr))');
    end
    X = cat(1, X, reshape(permute(Xi, [1 3 2]), npnts-p, nchs*p));
    Y = cat(1,Y,g.data(:, p+1:end,itr)');
end

% Indices for the diagonal elements (self connection)
diagIdx = vec(repmat((1:p*(nchs+1):p*nchs^2), p, 1) + repmat((0:(p-1))', 1, nchs))';
% Indices for the off-diagonal elements (connection to others)
offDiagIdx = setdiff(1:p*nchs^2, diagIdx);

% target vector
Y = vec(Y);
% predictor design matrix
X = blkdiageye(sparse(X),nchs);
X = [X(:, offDiagIdx) X(:, diagIdx)];
blks = [p*ones(1,nchs*(nchs-1)),p*nchs];

% normalize columns of X (requires CVX)
% if exist('norms','file')
%     X = X*spdiags(1./norms(X)',0,nchs^2*p,nchs^2*p); % K = sum(blks);
%     Y = Y./norm(Y);
% end


% if isempty(lambda)
%     % select lambda using heuristic
%     % Idea borrowed from Boyd et al [3]
%     if g.verb, fprintf('Using heuristic lambda selection...'); end
%         
%     cum_part = cumsum(blks(1:end-1));
%     
%     % guess regularization param for 
%     % each group of coefficients
%     K = length(blks)-1;
%     start_ind = 1;
%     lambdas = zeros(1,K);
%     
%     for i = 1:K,
%         sel = start_ind:cum_part(i);
%         lambdas(i) = norm(X(:,sel)'*Y);
%         start_ind = cum_part(i) + 1;
%     end
%     lambda_max = max(lambdas);
% 
%     % regularization parameter as fraction of 
%     % maximum group regularization parameter
%     g.admm_args.lambda = 0.1*lambda_max;   % 0.1*lambda_max
%     
%     if g.verb, fprintf('lambda set to %0.10g',lambda); end
% end


% Apply the ADMM method for group lasso estimation:
% group penalize AR coefficients with different time-lags
% between sources (nchs*(nchs-1) groups of size p). 
[initAR] = admm_gl(X, Y, blks, g.admm_args,'InitialState',initAR);  %,'InitialState',initAR

% assemble coefficient matrices
H2 = zeros(p,nchs,nchs);
H2(offDiagIdx) = vec(initAR(1:(p*nchs*(nchs-1))));          % non-diagonal elements
H2(diagIdx) = vec(initAR(((p*nchs*(nchs-1))+1):end));       % diagonal elements

AR = permute(full(H2), [3 2 1]);   % the recovered (nchs x nchs x p) connectivity matrix
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

