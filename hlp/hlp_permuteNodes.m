function CAT = hlp_permuteNodes(CAT,NodeOrder,allowDups)
% permute the order of nodes in a SIFT datastructure
%
% Author: Tim Mullen, SCCN/INC/UCSD 2014

if nargin<3
    allowDups = false;
end

if ~allowDups && length(unique(NodeOrder))~=length(NodeOrder)
    error('The elements of NodeOrder must be unique');
end

if length(NodeOrder)~=length(CAT.curComps)
    error('The length of NodeOrder must equal the number of nodes');
end

% permute the nodes
CAT = hlp_selectConnNodes(CAT,NodeOrder,[],[],false);

