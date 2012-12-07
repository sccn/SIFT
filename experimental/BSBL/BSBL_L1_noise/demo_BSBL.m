
clear; clc;

addpath('group_BP');

% problem dimension
M = 128;          % row number of the dictionary matrix 
N = 512;          % column number

blkNum = 7;       % active block number
blkLen = 8;        % block length

SNR = 15;
iterNum = 100;       % number of experiments

 

for it = 1: iterNum
    fprintf('\n\nRunning %d:\n',it);

    % Generate dictionary matrix with columns draw uniformly from the surface of a unit hypersphere
    Phi = randn(M,N);
    Phi = Phi./(ones(M,1)*sqrt(sum(Phi.^2)));
    
    % generate nonzero block coefficients
    beta = 1 - rand(blkNum,1)*0.2;     % correlation for each block
    %beta = ones(blkNum,1)*0.95;
    blks(:,1) = randn(blkNum,1);
    for i = 2 : blkLen
        blks(:,i) = beta .* blks(:,i-1) + sqrt(1-beta.^2).*(ones(blkNum,1).*randn(blkNum,1));
    end
    
    % normalize along row and vectorize
    nonzeroCoeff = reshape(blks',blkNum*blkLen,1);

    
    %========================================================================
    % put blocks at random locations and align with block partition (no overlap) 
    %========================================================================
    psbBlockNum = floor(N/blkLen);
    ind2 = randperm(psbBlockNum);
    indice2 = ind2(1:blkNum);
    Wgen = zeros(N,1);  blkLoc2 = [];
    for i = 1 : blkNum
        Wgen( (indice2(i)-1)*blkLen + 1 : indice2(i)*blkLen ) = nonzeroCoeff((i-1)*blkLen + 1 : i*blkLen);
        blkLoc2 = [blkLoc2,(indice2(i)-1)*blkLen + 1 : indice2(i)*blkLen];
    end
    
    % noiseless signal
    signal = Phi * Wgen;
    
    % Observation noise   
    stdnoise = std(signal)*10^(-SNR/20);
    noise = randn(M,1) * stdnoise;

    % Noisy signal
    Y = signal + noise;

 
    %======================================================================
    %              Benchmark (LS solution given support)
    %======================================================================
    supt = find( abs(Wgen)>0);
    x_ls = pinv(Phi(:,supt)) * Y;
    mse_bench(it) = (norm(Wgen(supt) - x_ls,'fro')/norm(Wgen,'fro'))^2; 
    fprintf('     bench :                             MSE: %g;\n', mean(mse_bench));
        
    
    %======================================================================
    %                       BSBL-L1
    %======================================================================
    blkStartLoc = [1:blkLen:N];
    iterNum = 3;          % [ Suggest setting: iterNum = 3 in general noisy cases; 
                          %                    iterNum = 6 in strongly noisy cases (SNR<5dB)]
    
    Result = BSBL_L1_noise(Phi,Y,blkLen,iterNum);
    
    for i = 1 : iterNum
        mse(it,i) = (norm(Wgen - Result.x(:,i),'fro')/norm(Wgen,'fro'))^2;
        time(it,i) = Result.time(i);
        fprintf('blockRewL1 : Iteration %d: Time: %5.4f;  MSE: %6.5f;\n', ...
            i,mean(time(:,i),1),mean(mse(:,i),1));
    end
    
    %======================================================================
    %                       BSBL-FM (0)
    %======================================================================
    learnLambda = 1;
    
    tic;
    Result6 = BSBL_FM(Phi,Y,blkStartLoc,learnLambda,'learnType',0);
    t_fm0(it) = toc;
    mse_fm0(it) = (norm(Wgen - Result6.x,'fro')/norm(Wgen,'fro'))^2;

    fprintf('\nBSBL-FM(ignore correlation) : time: %4.3f, MSE: %g',mean(t_fm0),mean(mse_fm0));
    
    
end



