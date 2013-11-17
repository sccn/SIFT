%% first, construct oscillators

clear

ndisc = 1000;   % samples to discard
% Nl = 100000+ndisc; %5000;
% Nr = 1; %100;
Nl = 5000+ndisc; %5000;
Nr = 100; %100;
q=3;
SRATE = 1;
f1=0.1; f2=0.4;     % oscillation freqs
c1=0.4; c2=0.8;    % coupling strength

a1 = 2*exp(-1/100)*cos(2*pi/10);
a2 = 2*exp(-1/100)*cos(2*pi/2.5);
b1 = -exp(-2/100);
b2 = b1;

A1 = [a1    0; 
      0     a2];        
A2 = [b1    0; 
      0     b2];
A=[A1 A2];
C = eye(q-1);
% clear A1 A2 A50;


%% Build AR model
q = size(C,1);
data = zeros(q,Nl,Nr);
for tr=1:Nr
   data(:,:,tr) = arsim(zeros(1,q),A,C,Nl)';
end
tmp = permute(data,[2 1 3]);    % time x chans x trials


% filter processes (single trial)
tmp(:,1) = eegfilt(data(1,:),1,f1-0.01,f1+0.01,0,30)';
tmp(:,2) = eegfilt(data(2,:),1,f2-0.01,f2+0.01,0,30)';

% normalize
tmp = zscore(tmp);
tmp(:,3) = zeros(size(tmp,1),1);
randvec=randn(size(tmp));

% construct third process
for t=50+1:Nl
    tmp(t,3) = c1*tmp(t-50,1)+c2*(tmp(t-50,2))+randvec(t-50);
end
tmp=tmp(ndisc+1:end,:,:);
    
%% Now create mixture model
data = zeros(size(tmp,1),2,size(tmp,3));
data(:,1,:) = sum(tmp(:,[1 2],:),2)+3*randn(size(tmp,1),size(tmp,3));
data(:,2,:) = tmp(:,3,:)+3*randn(size(tmp,1),size(tmp,3));

%% Now do the PDC