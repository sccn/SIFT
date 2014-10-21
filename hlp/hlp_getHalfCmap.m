function J = hlp_getHalfCmap(mapFcn,cmapSize,whichHalf)
% return the upper or lower half of a colormap
%
% Inputs:
%
%   mapFcn:    the name of the function that returns the colormap (e.g. 'jet'
%               or 'hsv'). Can also be a [N x 3] colormap matrix
%   cmapSize:   the cmapSize of the colormap (second argument to mapFcn function)
%   whichHalf:   'upper' or 'lower'
%
% Author: Tim Mullen, SCCN/INC/UCSD, May 2012

if nargin<2
    whichHalf = 'upper';
end

if isa(mapFcn,'function_handle') || ischar(mapFcn)
    J = feval(mapFcn,cmapSize);
else
    J = mapFcn;
end
switch whichHalf
    case 'upper'
        J = J(round(cmapSize/2):end,:);
    case 'lower'
        J = J(1:round(cmapSize/2)-1,:);
end
