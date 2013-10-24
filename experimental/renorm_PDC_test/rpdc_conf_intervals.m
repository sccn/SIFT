
% non-central confidence intervals (wierd plotting error)
% need to validate

pup = ncx2inv(0.95,2,size(data,1)*C.rPxy)/size(data,1);
pdn = max([zeros(size(pup)); 2*C.rPxy-pup]);
figure; errorbar(FOI,C.rPxy,pdn,pup);


pup = ncx2inv(0.95,2,size(data,1)*C.rPyx)/size(data,1);
pdn = max([zeros(size(pup)); 2*C.rPyx-pup]);
figure; errorbar(FOI,C.rPyx,pdn,pup);
