function [signal, state] = flt_loreta(varargin)
% Return the current source density for a given head model and data using
% the cortically-constrained standardized LORETA (low resolution electrical
% tomographic analysis) with a Bayesian update scheme for hyperparameters.
% The reconstructed CSD time-series (or source potential maps) will be 
% stored in signal.srcpot. This matrix has dimension [num_voxels x num_samples].
% 
% Author: Tim Mullen, Jan 2013, SCCN/INC/UCSD
%         Alejandro Ojeda, Jan 2013, SCCN/INC/UCSD
%         Christian Kothe, Jan 2013, SCCN/INC/UCSD

if ~exp_beginfun('filter'), return; end

declare_properties('name','LORETA','independent_channels',false, 'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg_nogui({'K','ForwardModel'},[],[],'Forward model (matrix)','shape','matrix'), ...
    arg_nogui({'L','LaplacianOperator'},[],[],'Laplacian operator. This is also known as the "precision matrix"'), ...
    arg_sub({'options','LoretaOptions'},{},...
    { ...
    arg({'maxTol','MaxTolerance'},1e-3,[0 Inf],'Tolerance for hyperparameter update loop','cat','Loreta Options'), ...
    arg({'maxIter','MaxIterations'},100,[1 Inf],'Maximum iterations for hyperparameter update loop','cat','Loreta Options'), ...
    arg({'verbose','VerboseOutput'},false,[],'Verbosity','cat','Loreta Options'), ...
    arg({'initNoiseFactor','InitialNoiseFactor'},0.001,[0 Inf],'Fraction of noise level. Used for initializing alpha parameter','cat','Loreta Options') ...
    arg({'block_size','BlockSize'}, [], [], 'Block granularity for processing. The inverse operator will be updated using blocks of this many samples. This assumes that the inverse solution is spatially stationary over this many samples.'), ...
    },'Additional options for Loreta function'), ...
    arg_nogui({'state','State'},[],[],'State object. When provided, hyperparameters will be estimated adaptively from prior state'));

if isempty(block_size)
    block_size = signal.pnts;
end

[C S] = size(signal.data);
numsplits = floor(S/block_size);

% convert options sub-struct to cell array
args = hlp_struct2varargin(options);

if isempty(state) || ~isfield(state,'s') || isempty(state.s)
    % mode is offline or we are initializing online filter
    % perform one-time SVD for faster computation.
    [U,S,V] = svd(K/L,'econ');
    state.iLV = L\V;
    state.s = diag(S);
    state.s2 = s.^2;
    state.Ut = U';
else
    % add current estimate of alpha, beta to args
    args = [args, state.alpha, state.beta];
end

% loop over sub-blocks and estimate CSD for each block
for i=0:numsplits-1
    range = 1+floor(i*S/numsplits) : min(S,floor((i+1)*S/numsplits));

    % call (bayesian) loreta estimator
    [signal.srcpot(:,range), state.alpha, state.beta, state.srcweights(:,:,i+1)] ...
        = dynamicLoreta(state.Ut, signal.data(:,range), state.s2, state.iLV, L, args{:});
end

% store the mean inverse operator            
state.srcweights = mean(state.srcweights,3);
    
exp_endfun;

