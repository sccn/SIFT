% addpath dal_ver1.05
% addpath prob
% addpath tprod

% samples
N = 500;

% channels
M = 20;

% trials
T = 10;

% MVAR model order
K = 3;

% number of interacting sources
nsources = 5;

% number of interactions
ninter = 3;

% dimensionality reduction
nPCA = 5;

% group lasso regularization parameter
lambda = 1e2;

% signal to noise ratio 
snr = 0.5;

% generate data
[S, H] = gen_ar_sech(nsources, (N+100)*T, K, ninter);
A = randn(M, nsources);
X = A*S;
X = reshape(X, M, N+100, T);
X = X(:, 1:N, :);
X = X / norm(vec(X));
noise = randn(size(X));
noise = noise / norm(vec(noise));
X1 = snr*X + (1-snr)*noise;

tic
[H_est1 PE1 out1] = mvar_scsa_em(struct('data', X1, 'morder', K, 'lambda', lambda, 'warmStart', [], 'PCA', nPCA));
toc

noise = randn(size(X));
noise = noise / norm(vec(noise));
X2 = snr*X + (1-snr)*noise;
tic
[H_est2 PE2 out2] = mvar_scsa_em(struct('data', X2, 'morder', K, 'lambda', lambda, ...
'warmStart', struct('B', out1.scsafilt, 'H', H_est1), 'PCA', struct('pcafilt', out1.pcafilt, 'pcapat', out1.pcapat)));
toc

