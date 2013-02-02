% The demo files are a mimicking of telemonitoring of fetal ECG recordings. 
% It re-produces the results in Section III.B of the following paper: 
%
% [1] Zhilin Zhang, Tzyy-Ping Jung, Scott Makeig, Bhaskar D. Rao, 
%     Compressed Sensing for Energy-Efficient Wireless Telemonitoring of 
%     Non-Invasive Fetal ECG via Block Sparse Bayesian Learning, submitted to 
%     IEEE Trans. on Biomedical Engineering. [Online] http://arxiv.org/abs/1205.1287
%
% Any question can be sent to Zhilin Zhang (zhangzlacademy@gmail.com)
%

clear;  close all;

% get the raw dataset (which includes strong noise)
load signal_01.mat;

% downsampling to 250 Hz
s = s(:,1:4:51200); 
 
% the size of the sensing matrix Phi
load Phi;
[M,N] = size(Phi);

% block size of the user-defined block partition of BSBL-BO
blkLen = 24;

% the user-defined block partition
groupStartLoc = 1:blkLen:N;

% variable for recovered dataset
X_hat = zeros(size(s));

%====================================================
% Compress all the ECG recordings, and recover them
%====================================================
k = 0;
for i = 1 : 8
    fprintf('\nChannel:%d\n',i);
    for j = 1 : size(s,2)/N
        k = k + 1;
        
        fprintf('  Segment %d: ',j);
        y = Phi * s(i,(j-1)*N+1:j*N)';   % sensing
        Y(i,(j-1)*M+1:j*M) = y';         % compressed data;
        
        % recover..
        tic;
        Result = BSBL_BO(Phi,y,groupStartLoc,0,'prune_gamma',-1, 'max_iters',20);
        runtime(k) = toc;
        
        X_hat(i,(j-1)*N+1:j*N) = ((Result.x))';
        
        mse(k) = (norm(s(i,(j-1)*N+1:j*N) - X_hat(i,(j-1)*N+1:j*N),'fro')/norm(s(i,(j-1)*N+1:j*N),'fro'))^2;
        fprintf(' MSE = %g, count = %d\n',mse(k),Result.count);

    end
    
end
fprintf('Total MSE: %g\n',mean(mse)); 

% display the original dataset 
set(0, 'DefaultFigurePosition', [100 70 500 600]);
ICAshow(s,'title','Original Recordings');

% display the reconstructed dataset 
set(0, 'DefaultFigurePosition', [650 70 500 600]);
ICAshow(X_hat,'title','Reconstructed Recordings by BSBL-BO');

save data_recovered_dataset;


% Next, we perform ICA decomposition of the recovered recordings, and
% perform the same ICA decomposition of the original recordings. Then, we
% compare the two ICA decompositions, seeing whether the extracted fetal
% ECG from the recovered recordings is the same as the extracted one from
% the original recordings.
%
% Note that in each running of the ICA algorithm, due to some randomness of  
% the algorithm, the fetal ECG could be the 3rd independent component, or
% the 4-th independent component
%

%====================================================
%   ICA decomposition of the original recordings
%====================================================

% band-pass filtering 
lo = 0.008;  % normalized low frequency cut-off
hi = 0.4;    % normalized high frequency cut-off
s_flt = BPFilter(s,lo,hi);

% remove mean and normalized to unit variance
Z = standarize(s_flt);

% whitening
[WhitenedSig,WhitenMatrix]=whiten(Z); 

% perform ICA decomposition (extracting 5 independent components)
s1 = FastICA(WhitenedSig,'method','defl','ICNum',5);
set(0, 'DefaultFigurePosition', [50 20 500 650]);

% display the ICA decomposition
seg = 1:1000;
ICAshow(s1(:,seg),'title','ICA of Original Recordings');



%====================================================
%   ICA decomposition of the recovered recordings
%====================================================

% perform the same band-pass filtering
X_hat_flt = BPFilter(X_hat,lo,hi);

% perform the same mean-removing and normalized to unit variance
Z_hat = standarize(X_hat_flt);

% perform the same whitening
[WhitenedSig2,WhitenMatrix2]=whiten(Z_hat); 

% perform the same ICA algorithm
[s2] = FastICA(WhitenedSig2,'method','defl','ICNum',5);


% sort the independent components such that they are in the same order of s1
s2_sort = zeros(size(s2));
for i = 1: size(s2,1)
    co = s1(i,:)*s2'/length(s1(i,:));
    [~,ind]=max(abs(co));
    s2_sort(i,:) = s2(ind,:)*sign(co(ind));
end
set(0, 'DefaultFigurePosition', [600 20 500 650]);

% Display the ICA decomposition from the recovered dataset.
% *** The fetal ECG could be the 3rd component, or the 4th component **** $
ICAshow(s2_sort(:,seg),'title','ICA of Recovered Recordings');

save data_final_results;


