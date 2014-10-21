function X = hlp_insertSingletonDim(X,dim)
% insert a singleton dimension into dimension 'dim' of matrix X
% author: Tim Mullen, 2011

nd = ndims(X);
v = [1:dim-1 nd+1 dim:nd];
[dummy, ind] = unique_bc(v,'first');
v = v(sort(ind));
X = permute(X,v);