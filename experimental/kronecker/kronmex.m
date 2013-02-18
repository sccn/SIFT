function c = kronmex(a,b)
%Kronecker product of the two matrix arguments in the order they appear
%If A is an m-by-n matrix and B is a p-by-q matrix, then the Kronecker
%product is the mp-by-nq block matrix
%   example, if X is 2 by 3, then KRONECKER(X,Y) is
%
%      [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
%        X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
%   Mayowa Aregbesola (Mayowa@ieee.org)
% Thanks to Etienne Grossmann (etienne@anonimo.isr.ist.utl.pt) & Bruno
% Luong <brunoluong@yahoo.com>

%Use the Mex function kronc
  c = kronc(a,b);
  