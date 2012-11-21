% Experiment demo for the first EEG experiment in: 
%   Zhilin Zhang, Tzyy-Ping Jung, Scott Makeig, Bhaskar D. Rao, 
%   Compressed Sensing of EEG for Wireless Telemonitoring with Low Energy 
%   Consumption and Inexpensive Hardware, accepted by IEEE Trans. on 
%   Biomedical Engineering, 2012
%
% This demo file shows the ability of BSBL-BO to recover signals which are
% non-sparse in the time domain or in any transformed domain
%
% Author: Zhilin Zhang (zhangzlacademy@gmail.com)
% Date  : Sep 12, 2012


clear; close all;

% load the first channel of the common EEG dataset in the EEGLab
load EEGdata_ch1;    

% load the randomly generated sparse binary matrix
% (you can use other sensing matrices)
load Phi;             
                     
% epoch length of the dataset is 384
N = 384;

% each epoch data is compressed to 192 data points (compression 50%)
M = 192; 

 
% We use DCT dictionary matrix
A=zeros(M,N);
for k=1:M
A(k,:)=dct(Phi(k,:));
end

recon_erp = []; 
count = 0;
for ch = 1 : 1       % channel number
    for ep = 1 : 80     % epoch number
        count = count + 1;
        
        x = EEGdata_ch1(ch,1:N,ep)';
        
        % compress an epoch
        y = Phi * double(x);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             use  BSBL-BO to recover
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        blkLen = 24;                  % the block size in the user-defined block partition
        groupStartLoc = 1:blkLen:N;   % user-defined block partition
        
        % run the algorithm
        tic;
        Result1 = BSBL_BO(A,y,groupStartLoc,0,'prune_gamma',-1, 'max_iters',7);
        t_bsbl(count) = toc;
        
        % restore recovered data
        recon_EEG(ch,1:N,ep) = (idct(Result1.x))';  
        
        % MSE
        mse_bsbl(count) = (norm(x - recon_EEG(ch,1:N,ep)','fro')/norm(x,'fro'))^2;
        
        % SSIM
        windowLen = 100;
        [mssim, ssim_map] = ssim_1d( x, recon_EEG(ch,1:N,ep)', windowLen);
        ssim_bsbl(count) = mssim;
        
        fprintf('BSBL-BO: time: %4.3f, MSE: %5.4f | SSIM: %5.4f\n',t_bsbl(count),mse_bsbl(count),ssim_bsbl(count));
        

    end
end


% plot the recovery result of an epoch EEG 
kepoch = 1;   
figure;
subplot(211); plot(EEGdata_ch1(1,:,kepoch),'linewidth',1); 
title('An Original EEG Epoch');
subplot(212); plot(recon_EEG(1,:,kepoch),'r','linewidth',1);
title_text = ['The Recovered EEG. MSE=',num2str(mse_bsbl(kepoch))];
title(title_text);

save recovery_result_by_BSBL;

