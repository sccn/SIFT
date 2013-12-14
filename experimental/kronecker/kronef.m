function c = kronef(a,b)
%Efficient Kronecker product of the two matrix arguments in the order they appear
%If A is an m-by-n matrix and B is a p-by-q matrix, then the Kronecker
%product is the mp-by-nq block matrix
%   example, if X is 2 by 3, then KRONECKER(X,Y) is
%
%      [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
%        X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
%   Mayowa Aregbesola
% Thanks to Etienne Grossmann (etienne@anonimo.isr.ist.utl.pt) & Bruno
% Luong <brunoluong@yahoo.com>


[ra, ca]=size(a);
  [rb, cb]=size(b);
  c = a(ones(rb,1)*(1:ra), ones(cb,1)*(1:ca)).* b((1:rb)'*ones(1,ra), (1:cb)'*ones(1,ca));
  