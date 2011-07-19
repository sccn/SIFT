function [AR PE ww] = mvar_dalSCSA(data,p,lambda,AR0,shrink_diagonal,loss,varargin)

% INPUTS:
%   data:       the data (nchs x npnts)
%   p:          the model order
%   lambda:     regularization coefficient
%   AR0:        initial solution (default: zeros)
%   varargin:   'name',value paired options to pass to dal.m
% OUTPUTS:
%   AR:         [nchs x nchs*p] VAR coefficient matrix
%   PE:         [nchs x nchs] estimated noise covariance matrix

% Modified by Tim Mullen, 2011 from s_test_hsgl.m from DAL 1.05
%
% This program infers a sparsely connected multivariate AR
% model assuming that the sources are independent and super
% Gaussian. We use the hyperbolic secant likelihood.
%
% In our IEEE TBME paper, we also estimate mixing matrix
% through an EM algorithm. This is omitted here for the sake
% of simplicity. Thus, all state variables are directly measured.
%
% Reference:
% Modeling sparse connectivity between underlying brain sources
% for EEG/MEG. Stefan Haufe, Ryota Tomioka, Guido Nolte,
% Klaus-Robert Mueller, and Motoaki Kawanabe, IEEE
% Trans. Biomed. Eng. 57(8), pp. 1954-1963, 2010.
% 
% Copyright(c) 2009-2011 Ryota Tomioka
%              2009      Stefan Haufe
% This software is distributed under the MIT license. See license.txt

persistent initAR;

[nchs npnts] = size(data);

if nargin<6
    loss = 'hs';
end
if nargin<5
    shrink_diagonal = true;
end

% initial solution
if nargin>4 && ~isempty(AR0)
    initAR = AR0;
elseif isempty(initAR)
    initAR = zeros(p*nchs^2,1);
end

X = [];
for jj = 1:p
  X = cat(3, X, data(:, p+1-jj:end-jj)'); 
end
X = reshape(permute(X, [1 3 2]), npnts-p, nchs*p);
% A = kron(speye(nchs), X);
A = blkdiageye(sparse(X),nchs);

%% Indices for the diagonal elements
indsB = vec(repmat((1:p*(nchs+1):p*nchs^2), p, 1) + repmat((0:(p-1))', 1, nchs))';

%% Indices for the off-diagonal elements
indsA = setdiff(1:p*nchs^2, indsB);

%% Design corresponding to the diagonal elements (self connection)
B = A(:, indsB);
if shrink_diagonal
    blks = [p*ones(1,nchs*(nchs-1)),p*nchs];
    Bu = [];
else
    Bu = B;
    B = [];
    blks = p*ones(1,nchs*(nchs-1));
end

%% Design corresponding to the off-diagonal elements (connection to others)
A = A(:, indsA);

%% Target
Y = vec(data(:, p+1:end)');

% Group penalize AR coefficients with different time-lags
% between sources (nchs*(nchs-1) groups of size p).
% Self connections are put into one group (size p*nchs).
if shrink_diagonal
    bias = []; %zeros(size(initAR));
    
    switch lower(loss)
        case 'sq'
            initAR=dalsqgl(initAR, [A B], Y, lambda,  'blks', blks,'display',1,varargin{:});
        case 'hs'
            [initAR] = dalhsgl(initAR,bias, [A B],[], ...
                               Y, lambda, 'blks', blks,'display',0,varargin{:});  % 'aa',zeros(size(A,1),1),
    end
else
   [initAR(indsA) initAR(indsB)] = dalhsgl(initAR(indsA),initAR(indsB), [A B],Bu, ...
        Y, lambda, 'blks', blks,'display',2,varargin{:});
end

H2 = zeros(p,nchs,nchs);
H2(indsA) = vec(initAR(1:(p*nchs*(nchs-1))));       % non-diagonal elements
H2(indsB) = vec(initAR(((p*nchs*(nchs-1))+1):end)); % diagonal elements

AR = permute(full(H2), [3 2 1]);   % the recovered (nchs x nchs x p) connectivity matrix
AR = reshape(AR,[nchs nchs*p]);

%% estimate noise covariance matrix
if nargout>1
    PE = cov(est_mvarResiduals(data,AR,zeros(1,nchs))');
end

if nargout>2
    ww = initAR;
end

function v = vec(x)
    v = x(:);
    
    
function C = blkdiageye(X,k)
% construct blockdiagonal matrix with k copies of X on diagonal
% Equivalent to C = kron(eye(k),X) but much faster

ss = repmat('X,',1,k); ss(end)=[]; 
eval(sprintf('C = blkdiag(%s);',ss));