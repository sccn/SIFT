%% first, construct oscillators

clear

ndisc = 1000;   % samples to discard
% Nl = 100000+ndisc; %5000;
% Nr = 1; %100;
Nl = 10000; %5000;
Nr = 1; %100;
q=3;
lag = [0 0 0;
       0 0 0;
       10 10 0];
       
SRATE = 128;
f1=10; f2=40;     % oscillation freqs
% c1=0.4; c2=0.8;    % coupling strength

tau = 100;           % damping time (samples)
beta = [0.04 0.04];  % smaller beta = smoother transition
c1max = 0.8;
c2max = 0.8;
c1 = c1max./(1+exp(-beta(1).*((1:Nl)-Nl/2)));
c2 = c2max./(1+exp(beta(2).*((1:Nl)-Nl/2)));

figure; plot(1:Nl,[c1;c2]);
legend('coupling 1->3','coupling 2->3');
vline(Nl,'--r');


a1 = 2*exp(-1/tau)*cos(2*pi*f1/SRATE);
a2 = 2*exp(-1/tau)*cos(2*pi*f2/SRATE);
b1 = -exp(-2/tau);
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
   data(:,:,tr) = arsim(zeros(1,q),A,C,[Nl Nr],ndisc)';
end
tmp = permute(data,[2 1 3]);    % time x chans x trials

%%
% filter processes (single trial)
tmp(:,1) = eegfilt(data(1,:),SRATE,f1-0.1,f1+0.1,0,30)';  
tmp(:,2) = eegfilt(data(2,:),SRATE,f2-0.1,f2+0.1,0,30)';

% normalize
tmp = zscore(tmp);

% initialize third process
tmp(:,3) = zeros(size(tmp,1),1);
randvec=randn(size(tmp));

% construct third process
% x3(t) = c1(t)x1(t-k1) + c2(t)x2(t-k2)
tmp(1:max(lag(:)),3) = 0;

for t=max(lag(:))+1:Nl
    tmp(t,3) = c1(t)*tmp(t-lag(3,1),1)+c2(t)*(tmp(t-lag(3,2),2));  %+randvec(t-lag(3,2))
end
% tmp=tmp(ndisc+1:end,:,:);
    
eegplot(tmp','srate',SRATE);

%% Now create mixture model
noise_std = 1; %3
data = zeros(size(tmp,1),2,size(tmp,3));
data(:,1,:) = sum(tmp(:,[1 2],:),2)+noise_std*randn(size(tmp,1),size(tmp,3));   % y1 = x1 + x2 + noise
data(:,2,:) = tmp(:,3,:)+noise_std*randn(size(tmp,1),size(tmp,3));              % y2 = x3 + noise

eegplot(data','srate',SRATE);

%% plot the spectrum
figure; 
subplot(211); spectrogram(data(:,1),1000,900,4096,SRATE);
subplot(212); spectrogram(data(:,2),1000,900,4096,SRATE);

%% Now do the PDC