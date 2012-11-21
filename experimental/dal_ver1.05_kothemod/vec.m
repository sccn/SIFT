% vec - vectorize an array
%
% Copyright(c) 2009-2011 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function V=vec(M)
% V=reshape(M, [numel(M), 1]);
V = M(:);
