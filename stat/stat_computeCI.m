function ci = stat_computeCI(PConn,alpha,tail)
% compute confidence intervals across last dimension of matrix PConn
%
% Inputs:
% 
%       PConn:     connectivity distribution object returned from stat_bootstrap()
%                  or a matrix of distribution values for a single
%                  estimator
%       alpha:     alpha-significance threshold (e.g. alpha=0.05)
%       tail:      'upper': compute only upper ci (lower is set to mean of
%                           distribution)
%                  'lower': compute only lower ci (upper is set to mean of
%                           distribution)
%                  'both': compute upper and lower confidence intervals
%                          (alpha is divided by 2)
% Outputs:
%
%       ci:     the mean of the distribution
%
% See Also: stat_bootstrap()
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.8
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT
%
% Author: Tim Mullen, 2011, SCCN/INC, UCSD. 
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

if nargin<3
    error('SIFT:stat_computeCI','Insufficient number of arguments');
end

if isstruct(PConn)
    % multiple connectivity methods
    connmethods = hlp_getConnMethodNames(PConn);
    for m=1:length(connmethods)
        % recursively compute stats for all methods
        ci.(connmethods{m}) = stat_computeCI(PConn(cnd).(connmethods{m}),alpha,tail);  
    end
else
    % PConn is a matrix 
    sz = size(PConn);
    nd = length(sz);
    alpha = 100*alpha;
    
    switch lower(tail)
        case 'both'
            ci(1,:,:,:,:,:,:,:,:) = prctile(PConn,alpha/2,nd);      % lower
            ci(2,:,:,:,:,:,:,:,:) = prctile(PConn,100-alpha/2,nd);  % upper
        case 'upper'
            mval = mean(PConn,nd); % mean of estimator
            ci(1,:,:,:,:,:,:,:,:) = mval;
            ci(2,:,:,:,:,:,:,:,:) = prctile(PConn,100-alpha,nd);    % upper
        case 'lower'
            mval = mean(PConn,nd); % mean of estimator
            ci(1,:,:,:,:,:,:,:,:) = prctile(PConn,alpha,nd);        % lower
            ci(2,:,:,:,:,:,:,:,:) = mval;
        otherwise
            error('SIFT:stat_computeCI','unknown tail option');
    end
end
