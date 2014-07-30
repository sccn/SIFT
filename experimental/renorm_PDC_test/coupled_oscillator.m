%% first, construct oscillators

Nl = 100000; %5000;
Nr = 1; %100;
q=3;
SRATE = 1;
f1=0.1; f2=0.4;     % oscillation freqs
c1=0.4; c2=0.8;    % coupling strength

a1 = 2*exp(-1/100)*cos(2*pi/10);
a2 = 2*exp(-1/100)*cos(2*pi/2.5);
b1 = -exp(-2/100);
b2 = b1;

A1 = [a1    0   0; 
      0     a2  0
      0     0   0];        
A2 = [b1    0   0; 
      0     b2  0
      0     0  0];
A50 =[0     0   0
      0     0   0
      c1   c2   0];
A=[A1 A2];
A=[A1 A2 zeros(q,q*47) A50];
C = eye(q);
% clear A1 A2 A50;


%% Build AR model
q = size(C,1);
data = zeros(q,Nl,Nr);
for tr=1:Nr
   data(:,:,tr) = arsim(zeros(1,q),A,C,Nl)';
end
tmp = permute(data,[2 1 3]);    % time x chans x trials

%% Now create mixture model
data = zeros(size(tmp,1),2,size(tmp,3));
data(:,1,:) = sum(tmp(:,[1 2],:),2)+3*randn(size(tmp,1),size(tmp,3));
data(:,2,:) = tmp(:,3,:)+3*randn(size(tmp,1),size(tmp,3));

%% Now do the PDC