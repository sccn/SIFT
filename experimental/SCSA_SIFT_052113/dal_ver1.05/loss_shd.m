function varargout=loss_shd(zz, bb)
% max(abs(zz))

% zz = zz./max(abs(zz));

% zz(zz > 0.999) = 0.999;
% zz(zz < -0.999) = -0.999; 

gloss = (bb+atanh(zz));
floss = sum(zz.*(bb+atanh(zz))+log(sqrt(1-zz.^2)/pi));

% flossd = floss

hmin = 1;

if nargout<=3
  varargout={floss, gloss, hmin};
else
  hloss = spdiag(1./(1-zz.^2));
  varargout={floss, gloss, hloss, hmin};
%   if min(eig(hloss)) < 0
%     keyboard
%   end
end

% keyboard