% Example showing the ability of BSBL to recover non-sparse signals.
%
% The signal to recover is a real-world fetal ECG data, which consists two
% peaks of fetal ECG and one peak of maternal ECG.
%
% The goal is to recover the signal without distorting the two peaks of fetal
% ECG (the two peaks locate around at 65 and 180).
% 
% Details can be found in the paper:
%    Zhilin Zhang, Tzyy-Ping Jung, Scott Makeig, Bhaskar D. Rao, 
%    Compressed Sensing for Energy-Efficient Wireless Telemonitoring of 
%    Non-Invasive Fetal ECG via Block Sparse Bayesian Learning, submitted to 
%    IEEE Trans. on Biomedical Engineering, 2012. [Online] http://arxiv.org/abs/1205.1287
%
% Author: Zhilin Zhang (z4zhang@ucsd.edu)
% March 2012

clear; close all;  


N = 250;
M = 125; 
 

% load fetal ECG data
load ECGsegment.mat;
x = double(ecg);

% load a sparse binary matrix. The matrix was randomly generated, each
% column of which has only 15 entries of 1. Note that other similar sparse 
% binary matrix works well also.
load Phi.mat;

% =========================== Compress the ECG data ====================
y = Phi * x;


% ========================== Method 1 ===================================== 
%  Reconstruct the ECG data by BSBL-BO directly  
%==========================================================================

% Set block size. Other block sizes ranging from 10 to 30 also work well.
% You can try other values, e.g. 10, 15, 20, etc.
blkLen = 25;  

% Define block partition
groupStartLoc = 1:blkLen:N;
 
% Run BSBL-BO which exploits the intra-block correlation
% Note: The ECG data is non-sparse. To reconstruct non-sparse signals, you
% need to set the input argument 'prune_gamma' to any non-positive values
% (0, -1, -10, etc), namely turning off this pruning mechanism in SBL. But
% when you reconstruct sparse signals, you can use the default value (1e-2).
tic;
Result1 = BSBL_BO(Phi,y,groupStartLoc,0,'prune_gamma',-1, 'max_iters',20);
t_bo = toc;

set(0, 'DefaultFigurePosition', [20 150 500 400]);
figure(1);
subplot(211);plot(Result1.x); title('\bf Reconstructed by BSBL-BO Directly'); 
axis([-inf, inf, -50, 110]);
subplot(212);plot(x,'r');title('\bf Original ECG Signal (Peaks at 65 and 180 are fetal ECG)'); 
axis([-inf, inf, -50, 110]);

disp(['Direct Method :',num2str(10*log10((norm(x-Result1.x)^2)/norm(x)^2))]);

% =========================== Second Method ==============================
% First recover the signal's coefficients in the DCT domain;
% Then recover the signal using the DCT ceofficients and the DCT basis
%=========================================================================

% Look at the coefficients when representing the fetal ECG signal in DCT
% domain; Still, the coefficients are not sparse. To recover all the
% coefficients are not impossible for most compressed sensing algorithms!
set(0, 'DefaultFigurePosition', [400 150 500 400]);
figure(2);
plot(dct(ecg)); title('\bf DCT coefficients of the fetal ECG signal; They are still not sparse!');

% Now we use the BSBL-BO. Its' block size is randomly set ( = 25, as in the
% above experiment).
A=zeros(M,N);
for k=1:M
A(k,:)=dct(Phi(k,:));
end
tic;
Result2 = BSBL_BO(A,y,groupStartLoc,0,'prune_gamma',-1, 'max_iters',20);
t_bo2 = toc;

signal_hat = idct(Result2.x);

set(0, 'DefaultFigurePosition', [800 150 500 400]);
figure(3);
subplot(211);plot(signal_hat); title('\bf Reconstructed by BSBL-BO from DCT Domain'); 
axis([-inf, inf, -50, 110]);
subplot(212);plot(x,'r');title('\bf Original ECG Signal (Peaks at 65 and 180 are fetal ECG)'); 
axis([-inf, inf, -50, 110]);

disp(['With DCT      :',num2str(10*log10((norm(x-signal_hat)^2)/norm(x)^2))]);
