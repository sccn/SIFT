% This demo is to show the powerful ability of SBL algorithms at the situation when
% the sensing matrix is highly coherent.
%
% This demo mimics an EEG source localization scenario, where the sensing matrix is 
% a simplified lead-field matrix construncted from a brain model. The
% correlation among columns of the sensing matrix is more than 0.99.
%
% In most trials, T-MSBL is better than M-SBL, and much better than
% M-FOCUSS. In fact, T-MSBL is also much better than most other existing
% MMV algorithms.
%
% ***********************
% Please note that in the case when sensing matrices are highly coherent,
% the input argument settings are slightly different from those in other
% cases.
% ***********************
%
% Any questions about this demo can be sent to zhangzlacademy@gmail.com
%
%


clear; clc;  close all;

L = 5;           % snapshot number
M = 400;         % column number of the dictionary matrix
N = 50;          % row number of the dictionary matrix
K = 3;           % number of the nonzero coefficients at the initial stage



load LF;   % load the sensing matrix (the so-called lead matrix constructed from a brain model)


% only use a simplified version of the sensing matrix
Phi = LF(1:N,1:M);
Phi = Phi./(ones(N,1)*sqrt(sum(Phi.^2)));

% Generate the K nonzero rows, each row being an AR(1) process. All the AR(1)
% processes have different AR coefficients, which are randomly chosen from [0.8,1)
beta = 1- rand(K,1)*0.2;           

nonzeroW(:,1) = randn(K,1);
for i = 2 : L*100
    nonzeroW(:,i) = beta .* nonzeroW(:,i-1) + sqrt(1-beta.^2).*(ones(K,1).*randn(K,1));
end
nonzeroW = nonzeroW(:,end-L+1:end);   % Ensure the AR processes are stable


% Rescale each row such that the squared row-norm distributes in [1,scalefactor]
nonzeroW = nonzeroW./( sqrt(sum(nonzeroW.^2,2)) * ones(1,L) );
scalefactor = 3;
mag = rand(1,K); mag = mag - min(mag);
mag = mag/(max(mag))*(scalefactor-1) + 1;
nonzeroW = diag(sqrt(mag)) * nonzeroW;

% Locations of nonzero rows are randomly chosen
ind = randperm(M);
indice = ind(1:K);
Wgen = zeros(M,L);
Wgen(indice,:) = nonzeroW;

% Noiseless signal
signal = Phi * Wgen;


% observation noise
stdnoise = 0.02;
noise = randn(N,L) * stdnoise;

% noisy signal
Y = signal + noise;

% compute total SNR
SNR_rec = 20*log10(norm(signal,'fro')/norm(noise,'fro'));
fprintf('\nSNR = %g dB  \n',SNR_rec);


 
%==================================================
%      support
%==================================================
supp = find(mean(abs(Wgen),2)>0);


%============================================
%             T-MSBL
%============================================
b = -0.9;
for kk = 1 : L, bvect(kk) = b^(kk-1); end
B0 = toeplitz(bvect);

% Note that in the situation when sensing matrices are highly coherent, 
% the input argument settings are slightly different from those in other cases. 
% You need to fix B to a fixed value, such as B0 here. In a lot of experiments 
% on real-world data when sensing matrices are highly coherent, we also found 
% this setting of B0 was helpful. 
[X_sbl3] = TMSBL(Phi, Y, 'fix_B',B0,'enhance_lambda',1);

F_sbl3 = perfSupp(X_sbl3,supp,'firstlargest', K);
Finx_sbl = (F_sbl3==1);


err_total_sbl3 = (norm(Wgen - X_sbl3,'fro')/norm(Wgen,'fro'))^2;
fprintf('T-MSBL(unknown noise variance): MSE = %g | Success Rate = %g(=1: perfect)\n',err_total_sbl3,Finx_sbl);



% ===========================================
%          M-SBL (given the true noise variance)
%============================================

[X_m] = MSBL(Phi, Y, stdnoise^2, 0);
F_m = perfSupp(X_m,supp,'firstlargest', K);
Finx_m = (F_m==1);

err_total_m = (norm(Wgen - X_m,'fro')/norm(Wgen,'fro'))^2;
fprintf('M-SBL ( given true noise var.): MSE = %g | Success Rate = %g(=1: perfect)\n',err_total_m,Finx_m);


% ===========================================
%          MFOCUSS (given the true noise variance)
%============================================

[X_tf] = MFOCUSS(Phi, Y, stdnoise^2);
F_tf = perfSupp(X_tf,supp,'firstlargest', K);
Finx_tf = (F_tf==1);


err_total_tf = (norm(Wgen - X_tf,'fro')/norm(Wgen,'fro'))^2;
fprintf('MFOCUSS(given true noise var.): MSE = %g | Success Rate = %g(=1: perfect)\n',err_total_tf,Finx_tf);





