function res = hlp_findscaling(X, scaling)
% Obtain information necessary to scale the given data. Works with hlp_applyscaling.
% Scale-Info = hlp_findscaling(Data, Scale-Mode)
%
% This is just a convenience tool to implement simple data (e.g. feature) scaling operations.
%
% In:
%   Data        : data matrix of [Observations x Variables]
%   Scale-Mode  : scaling mode, one of {std,minmax,whiten}
%
% Out:
%   Scale-Info  : scaling structure that can be used with hlp_applyscaling, to scale data
%
%               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%               2010-03-28

if ~exist('scaling','var') 
    scaling = 'minmax'; end

switch scaling
    case 'center'
        res = struct('add',{-mean(X)});
    case 'std'
        res = struct('add',{-mean(X)}, 'mul',{1./ std(X)});
    case 'minmax'
        res = struct('add',{-min(X)}, 'mul',{1./ (max(X) - min(X))});
    case 'whiten'
        [Uc,Lc] = eig(cov(X));
        res = struct('add',{-mean(X)},'project',{Uc * sqrt(inv(Lc))'});
    otherwise
        if ~isempty(scaling) && ~strcmp(scaling,'none')
            error('hlp_findscaling: unknown scaling mode specified'); end
        res = struct();
end