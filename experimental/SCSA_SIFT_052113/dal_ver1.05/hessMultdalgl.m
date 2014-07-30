% hessMultdalgl - function that computes H*x for DAL with grouped
%                 L1 regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function yy = hessMultdalgl(xx, A, eta, Hinfo)

blks =Hinfo.blks;
hloss=Hinfo.hloss;
I    =Hinfo.I;
vv   =Hinfo.vv;
nm   =Hinfo.nm;
lambda=Hinfo.lambda;

yy = hloss*xx;

nn=sum(blks);

cumblks=[0,cumsum(blks)];
for kk=1:length(I)
  jj=I(kk);
  bsz=blks(jj);
  J=cumblks(jj)+(1:bsz);
  vn=vv(J)/nm(jj);

  ff=lambda/nm(jj);
  AJ=A.slice(J);
  xk=AJ'*xx;
  yy = yy + eta(1)*(AJ*((1-ff)*xk+ff*(vn'*xk)*vn));
end

B=Hinfo.B;
if ~isempty(B)
  yy = yy + eta(2)*(B*(B'*xx));
end
