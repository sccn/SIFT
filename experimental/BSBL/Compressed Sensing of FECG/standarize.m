function standSig = standarize(OldSig)
% standarize --- standarize the input signals, i.e., make each signal zero
% mean and unit variance. If the input OldSig is in the form N*P, where N 
% is the number of signals, and P is the length of each signal, then returned 
% standSig is also in the form N*P, and each row (each signal) has zero
% mean and unit variance.
%
% Command:
%     standSig = standarize(OldSig)
%
% See also:
%     remstd    centering   whiten    
%
% Author:Zhi-Lin Zhang
%
% version: 1.0     Date: Apr.17,2005



sig = OldSig - mean (OldSig')' * ones (1,size (OldSig, 2));

for t = 1: size(sig,1)
    standSig(t,:) = sig(t,:)/std(sig(t,:));
end



