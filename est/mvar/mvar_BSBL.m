function [AR PE State] = mvar_BSBL(varargin)
% Algorithm: BSBL
%
% Description:
%
% BSBL Algorithm by Zhilin Zhang
%
% References and Code:
%
% Dependencies: BSBL_BO()
%
% ------------------------------------------------------------------------
% INPUTS:
%   data:       the data (nchs x npnts)
%   p:          the model order
%   lambda:     regularization parameter
%   rho:        the augmented Lagrangian parameter
%   alpha:      the over-relaxation parameter (typical 
%               values for alpha are between 1.0 and 1.8)
%   initState:        initial solution (default: zeros)
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
                arg({'normcols','NormCols'},'zscore',{'none','norm','zscore'},'Normalize columns of dictionary'), ...
                arg_subtoggle({'warmStart','WarmStart'},[], ...
                {...
                    arg({'initState','InitialState'},[],[],'Initial BSBL state object. If empty, defaults will be used.') ...
                },'Warm start. The previously estimated state will be used as a starting estimate for successive BSBL operations. Alternately, you may provide an non-empty initial state structure via the ''InitialState'' argument.'), ...
                arg_sub({'bsbl_args','BSBL_Options'},[],{ ...
                       arg({'LearnLambda','LambdaLearningRule'},'LowNoise',{'off','LowNoise','HighNoise'},{'Lambda learning rule',sprintf(['\n' ...
                                                                                                              '''LowNoise'' :  Use the lambda learning rule for very LOW SNR cases (SNR<10dB). Uses lambda=std(y)*1e-2 or user-input value as initial value.\n\n' ...
                                                                                                              '''HighNoise'':  Use the lambda learning rule for medium noisy cases (SNR>10dB). Uses lambda=std(y)*1e-2 or user-input value as initial value.\n\n' ...
                                                                                                              '''off'':        Do not use the lambda learning rule. Uses lambda=1e-14 or user-input value as initial value' ...
                                                                                                              ])}) ...
                       arg({'PRUNE_GAMMA','GammaPruningThreshold'},[],[0 Inf],{'Gamma Pruning threshold',sprintf(['\n' ...
                                                                                   'Threshold for prunning small hyperparameters gamma_i. Blocks of parameters are pruned when their ''power'' (gamma_i) is small.\n'    ...
                                                                                   'In noisy cases, you can set PRUNE_GAMMA = 1e-3 or 1e-4. \n'   ...
                                                                                   'In strong noisy cases (e.g. SNR <= 6 dB), set PRUNE_GAMMA = 1e-2 for better performance.' ...
                                                                                 ])}), ...
                       arg({'LAMBDA','InitialLambda','Lambda'},[],[],'Initial regularization value (lambda). If empty, the initial value is estimated from the data.','type','denserealdouble'), ...                                                        
                       arg({'EPSILON','StoppingTolerance'},1e-8,[0 Inf],'Stopping tolerance. The algorithm terminates when the L2 norm of the change in parameter estimates is small than this value.'), ...
                       arg({'LEARNTYPE','LearnCorrelation'},true,[],'Exploit correlation structure.'), ...
                       arg({'MAX_ITERS','MaxIterations'},600,[1 Inf],'Maximum iterations'), ...
                       arg({'PRINT','VerboseOutput'},false,[],'Verbosity') ...
                   },'Additional options for BSBL algorithm') ...
            );
                
% arg_toworkspace(g);
  
% modify some arguments
switch g.bsbl_args.LearnLambda
    case 'off',         g.bsbl_args.LearnLambda = 0;
    case 'LowNoise',    g.bsbl_args.LearnLambda = 1;
    case 'HighNoise',   g.bsbl_args.LearnLambda = 2;
end

[nchs npnts ntr] = size(g.data);
p = g.morder;

% initialize state
if g.warmStart.arg_selection && ~isempty(g.warmStart.initState)
    initAR = g.warmStart.initState;
elseif ~g.warmStart.arg_selection
    % reset initAR
    initAR = struct([]);
end
if isfield(initAR,'x') && size(initAR.x,1) ~= p*nchs^2
    % dimensions have changed, reset state
    fprintf('mvar_BSBL: model dimensions changed -- resetting state\n');
    initAR = struct([]);
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
X = [X(:, offDiagIdx) X(:, diagIdx)];
blks = [p*ones(1,nchs*(nchs-1)),p*nchs];
blks = [1 cumsum(blks(1:end-1))+1];



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


% Apply the BSBL method for sparse VAR estimation:
% group penalize AR coefficients with different time-lags
% between sources (nchs*(nchs-1) groups of size p). 
g.bsbl_args.InitState = initAR;
args = hlp_struct2varargin(g.bsbl_args,'suppress',{'arg_direct','LearnLambda'});
[initAR] = BSBL_BO(X,Y,blks,g.bsbl_args.LearnLambda,args{:});

% keyboard;
% [initAR] = admm_gl(X, Y, blks, g.admm_args,'InitialState',initAR);  %,'InitialState',initAR

% assemble coefficient matrices
H2 = zeros(p,nchs,nchs);
H2(offDiagIdx) = vec(initAR.x(1:(p*nchs*(nchs-1))));          % non-diagonal elements
H2(diagIdx) = vec(initAR.x(((p*nchs*(nchs-1))+1):end));       % diagonal elements

AR = permute(full(H2), [3 2 1]);   % the recovered (nchs x nchs x p) connectivity matrix
AR = reshape(AR,[nchs nchs*p]);

%% estimate noise covariance matrix
if nargout>1
    res = est_mvarResiduals(g.data,AR,zeros(1,nchs));
    res = res(:,:);
    PE = cov(res');
end

if nargout>2
    State = initAR;
end

function v = vec(x)
v = x(:);

