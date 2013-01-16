function [AR PE lambdaOpt] = mvar_ridge(varargin)
% Algorithm: Ridge Regression
%
% Description:
%
% This method infers a fully connected multivariate
% autoregressive (VAR) model using L2 regularization,
% otherwise known as "ridge-regression", "Tikhonov 
% regularization" or "minimum-L2-norm" [2]
%
% VAR[p] coefficients inferred using ridge regression
% (L2 regularization). We assume most VAR coefficients
% (interactions) are small, and apply an L2 penalty
% for large coefficients. A consequence is that VAR
% coefficients are never exactly zero; thus statistical
% thresholding is an important follow-up step.
% This constraint allows us to solve under-determined
% systems (e.g. more parameters than observations).
%
% If Y = D(y,p) is a p-lag delay embedding of multi-
% variate data vector y(t), X is a block toeplitz 
% matrix of lagged copies of the delay embedded data,
% and A is an augmented matrix of VAR coefficients,
% then we may adopt the structural model:
%
% Y = XA + e,  for gaussian noise e
%
% We then seek to solve the optimization problem
%
% A_hat = argmin_A{0.5||Y-XA||_2^2 + L*||A||_2}
%
% where L is the regularization parameter. 
% L can be determined automatically using Generalized 
% Cross Validation (GCV) [3].
%
% This algorithm is a good choice if you have few 
% data samples and a fairly large number of sources
% and/or a high model order. Note that, unlike sparse
% methods, ridge regression yields a fully connected
% graph (all VAR coefficients are non-zero). Follow up
% with statistical thresholding.
%
% Dependencies: ridgeGCV()
%
% ------------------------------------------------------------------------
% INPUTS:
%   data:       the data (nchs x npnts)
%   p:          the model order
%   lambda:     regularization coefficient
% OUTPUTS:
%   AR:         [nchs x nchs*p] VAR coefficient matrix
%   PE:         [nchs x nchs] estimated noise covariance matrix
%
% See Also: est_fitMVAR(), mvar_dalSCSA(), est_fitMVARKalman()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual.
%     Available at: http://www.sccn.ucsd.edu/wiki/SIFT
% [2] Hoerl, A.E.; R.W. Kennard (1970). Ridge regression: Biased estimation 
%     for nonorthogonal problems. Technometrics 42(1):80?86. JSTOR 1271436.
% [3] Golub G, Heath M, Wahba G. (1979) Generalized Cross-validation as a  
%     method for choosing a good ridge parameter. Technometrics. 
%     doi: 10.1080/00401706.1979.10489751
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

verb = arg_extract(varargin,{'verb','Verbosity'},[],0);

g = arg_define([0 1],varargin, ...
                arg_norep({'data','Data'},mandatory,[],'Data Matrix. Dimensions are [nchs x npnts x ntr].'), ...
                arg({'p','ModelOrder','morder'},10,[],'VAR Model order'), ...
                arg_subswitch({'lambdaSelMode','LambdaSelectionMode'},'automatic', ...
                {'manual',{ ...
                    arg({'lambda','RegularizationParam'},1,[0 Inf],'Regularization parameter (lambda)','type','denserealdouble') ...
                    }, ...
                 'automatic',{ ...
                    arg({'lambdaGridSize','LambdaGridSize'},100,[0 Inf],'Grid size for regularization param search. This is used to automatically select the regularization parameter which minimizes the Generalized Cross-Validation (GCV) criteria.') ...
                    } ...
                 },'Selection mode for lambda. Automatic (GCV grid search) or Manual (must provide lambda)'), ...
                arg({'varPrior','VariancePrior'},1,[0 Inf],'Parameter covariance prior (sigma). If sigma a scalar, the parameter covariance prior is taken to be a diagonal matrix with sigma (variance) on diagonals. Otherwise, sigma can be a full prior covariance matrix. A sparse matrix is advised if covariance matrix is not dense.'), ...
                arg({'normcols','NormCols'},'zscore',{'none','norm','zscore'},'Normalize data. Columns of data matrix X and target vector Y are normalized using the chosen method.'), ...
                arg_subtoggle({'warmStart','WarmStart'},[], ...
                {...
                    arg({'initState','InitialState'},[],[],'Initial BSBL state object. If empty, last state will be used.') ...
                },'Warm start. The previously estimated state will be used as a starting estimate for successive BSBL operations. Alternately, you may provide an non-empty initial state structure via the ''InitialState'' argument.'), ...
                arg({'verb','Verbosity'},verb,[],'Verbose output','type','int32','mapper',@(x)int32(x)) ...
                );

% arg({'scaled','UseScaling'},true,[],'Scaling option. If set, coefficient estimates are restored to the scale of the original data'), ...

arg_toworkspace(g);

[nchs npnts ntr] = size(data);

% % initialize state
% if g.warmStart.arg_selection && ~isempty(g.warmStart.initState)
%     initAR = g.warmStart.initState;
% elseif ~g.warmStart.arg_selection
%     % reset initAR
%     initAR = struct([]);
% end
% if isfield(initAR,'x') && size(initAR.x,1) ~= p*nchs^2
%     % dimensions have changed, reset state
%     fprintf('mvar_BSBL: model dimensions changed -- resetting state\n');
%     initAR = struct([]);
% end

% assemble the predictor (design) matrix and target vector
[X Y] = hlp_mkVarPredMatrix(g.data,p);

% normalization
switch lower(g.normcols)
    case 'zscore'
        % standardize columns of X and Y (unit variance)
        X = zscore(X,1,1);
        Y = zscore(Y,1,1);
    case 'norm'
        % normalize columns of X and Y (unit norm)
        X = bsxfun(@rdivide,X,sqrt(sum(X.^2)));
        Y = bsxfun(@rdivide,Y,sqrt(sum(Y.^2)));
end

% It is convenient to express A as a column vector so we must...
% vectorize the target matrix...
Y = hlp_vec(Y);
% ...and expand predictor matrix
X = blkdiageye(sparse(X),nchs);


switch g.lambdaSelMode.arg_selection
    case 'manual'
        lambdaOpt = g.lambdaSelMode.lambda;
        gridSize  = 0;
    case 'automatic'
        gridSize  = g.lambdaSelMode.lambdaGridSize;
        lambdaOpt = [];
end
       
if isscalar(g.varPrior)
    g.varPrior=g.varPrior*speye(size(X,2));
end

[H2 lambdaOpt] = ridgeGCV(Y,X,g.varPrior,lambdaOpt,gridSize,false,verb);

H2 = reshape(H2,[p,nchs,nchs]);

AR = permute(full(H2), [3 2 1]);   % the recovered (nchs x nchs x p) connectivity matrix
AR = reshape(AR,[nchs nchs*p]);

%% estimate noise covariance matrix
if nargout>1
    res = est_mvarResiduals(data,AR,zeros(1,nchs));
    res = res(:,:);
    PE = cov(res');
end
