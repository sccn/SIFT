function varargout = sim_simulateSources(varargin)
% Project source amplitudes through a forward model to generate channel data.
%
% This function is a wrapper for sim_fwdProj() and maintained for backwards-
% compatibility. This function may be deprecated in a future release.
% 
% Author: Tim Mullen, SCCN/INC/UCSD, 2013

[arg1 arg2 arg3 arg4 arg5 arg6 arg7] = sim_fwdProj(varargin{:});
for k=1:nargout
    eval(sprintf('varargout{%d} = arg%d;',k,k));
end