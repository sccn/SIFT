% gl_softth - soft threshold function for grouped L1 regularization
% 
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [vv,ss]=gl_softth(vv, lambda,info)

if all(info.blks==info.blks(1))
  n=length(vv);
  bsz=info.blks(1);
  vv=reshape(vv,[bsz,n/bsz]);
  ss0=sqrt(sum(vv.^2));
  ss=max(ss0-lambda,0);
  vv=vv*spdiag(ss./ss0);
  vv=vv(:);
else
  ss=zeros(length(info.blks),1);
  ix0=0;
  for kk=1:length(info.blks)
    I=ix0+(1:info.blks(kk));
    ix0=I(end);
    ss(kk)=norm(vv(I));
    ssk=max(ss(kk)-lambda,0);
    if ssk>0
      vv(I)=ssk/ss(kk)*vv(I);
    else
      vv(I)=0;
    end
    ss(kk)=ssk;
  end
end


