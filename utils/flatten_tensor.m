
function Bf = flatten_tensor(B,inds)
% Bf = flatten_tensor(B,inds)
% "flatten" a 1-D cell array of [M_i x M_i x Q] tensors into a
% 1-D cell array of [1 x Q] vectors
%
% Inputs:
%
% B is a 1:S cell array where
% B{si} = [M_i x M_i x Q] connectivity tensor.
% We convert connectivity to cell array of form:
% Bf = {s1(1,1) s1(1,2) s1(1,3) ... s1(N,1) s1(N,2) s1(N,3) ...
%       s2(1,1) s2(1,2) s2(1,3) ... }
% where si(i,j) is the 1 x Q vector of connectivity values from channels 
% j to i, for the ith subject 
% 
% inds{si} determines which connectivity pairs to retain in the flattened
% array. If this is omitted, then all values are retained.
% inds{si} is a list of linear indices into the first two dims of B{si}
% That is, inds{si}(k) indexes to some B{si}(i,j,:) where 
% k = sub2ind([M_i, M_i],i,j)
% Linear indices for partitions can be created using sub2ind(). 
% For instance, the following code separates the diagonal and off-diagonal
% linear indices of B{si}:
%
% for si=1:length(B)
%     % get linear indices of diagonal and off-diagonal elements
%     N = size(B{si},1);
%     all_inds = 1:N^2;
%     dg_inds{si} = sub2ind([N N],1:N,1:N);       % <-- diagonal idx
%     od_inds{si} = setdiff_bc(all_inds,dg_inds); % <-- off-diagonal idx
% end
% 
% Output:
%
% Bf is the flattened tensor, in the format noted above
%
% Author: Tim Mullen, 2014
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

if nargin < 2, inds = []; end
nsubj = length(B);
Bf = {};
for si=1:nsubj
    % make 2-D
    sz = size(B{si});
    B{si} = reshape(permute(B{si},[2 1 3]),prod(sz(1:2)),sz(3));
    % determine indices for this subject
    if isempty(inds), inds_si = 1:sz(1)^2;
    else inds_si = inds{si};
    end
    % select rows (note this is selecting elements of B{si}(i,j,:)
    B{si} = B{si}(inds_si,:);
    % convert each row to an element of cell array
    Bf = [Bf, mat2cell(B{si},ones(1,size(B{si},1)),size(B{si},2))'];
end

