
function NodeValue = hlp_computeGraphMeasure(causality,ch1,selectedvars,graphMeasure)
% NodeValue = hlp_computeGraphMeasure(causality,ch1,selectedvars,graphMeasure)
%
% Compute a univariate graph-theoretic measure for a given node of a graph
%
% Inputs:
%
%   causality:      [num_vars x num_vars x <num_times>] causal matrix
%                   obtained from est_mvarConnectivity() collapsed across
%                   frequencies
%   ch1:            index of node
%   selectedvars:   indices of other nodes
%   graphMeasure:   which measure to compute (see below)
%
% Outputs:
%   
%   NodeValue:      the graph measure for the desired node (ch1)
%
%
% See Also: vis_causalBrainMovie3D(), est_mvarConnectivity(),
%           hlp_collapseFrequencies()
%
%
% References: 
% 
%   Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT
%
% Author: Tim Mullen, 2010, SCCN/INC, UCSD. 
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


% get indices of all nodes we will plot, except the current one
othervars = setdiff(selectedvars,ch1);

switch lower(graphMeasure)
    case 'none'
        NodeValue = zeros(1,size(causality,3));
    case 'outflow'
        % Compute outflow from ch1 in each freq
        NodeValue = squeeze(sum(causality(othervars,ch1,:),1));
    case 'mag_outflow'
        % Compute outflow from ch1 in each freq, ignoring sign
        NodeValue = squeeze(sum(abs(causality(othervars,ch1,:)),1));
    case 'inflow'
        % Compute inflow to ch1 in each freq
        NodeValue = squeeze(sum(causality(ch1,othervars,:),2));
    case 'causalflow'
        outflow =   squeeze(sum(causality(othervars,ch1,:),1));
        inflow  =   squeeze(sum(causality(ch1,othervars,:),2));
        NodeValue = outflow - inflow;
    case 'outdegree'
        % Compute number of outgoing edges from ch1 in each freq
        NodeValue = squeeze(sum(logical(causality(othervars,ch1,:)),1));
    case 'indegree'
        % number of incoming edges to ch1 in each freq
        NodeValue = squeeze(sum(logical(causality(ch1,othervars,:)),2));
    case 'causaldegree'
        % outdegree - indegree 
        outflow =   squeeze(sum(logical(causality(othervars,ch1,:)),1));
        inflow  =   squeeze(sum(logical(causality(ch1,othervars,:)),2));
        NodeValue = outflow - inflow;
    case 'asymmetryratio'
        % 1 if all edges are outgoing, -1 if all edges are incoming.
        % 0 if balanced
        outflow =   squeeze(sum(causality(othervars,ch1,:),1));
        inflow  =   squeeze(sum(causality(ch1,othervars,:),2));
        NodeValue = (outflow - inflow)./(outflow+inflow);
    otherwise
        % user wants to map a different Conn measure to this
        % (e.g., ERSP)
end