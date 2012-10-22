function X = hlp_applyscaling(X, si)
% Scale the given data according to some info, as obtained from hlp_findscaling.
% Out-Data = hlp_applyscaling(In-Data, Scale-Info)
% 
% This is just a convenience tool to implement simple data (e.g. feature) scaling operations.
%
%                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-03-28

if isfield(si,'add') X = X+repmat(si.add,[size(X,1),1]); end
if isfield(si,'mul') X = X.*repmat(si.mul,[size(X,1),1]); end
if isfield(si,'project') X = X*si.project; end
if isfield(si,{'add','mul','project'}) X(~isfinite(X(:))) = 0; end