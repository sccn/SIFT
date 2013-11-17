% This demo shows how to set input parameters of BSBL-EM and BSBL-BO 
% for noisy experiments when block partition is known.
%
% Also it shows the performance improvement if exploiting intra-block correlation. 


clear; 

% problem dimension
M = 80;          % row number of the dictionary matrix 
N = 162;          % column number

blkNum = 5;       % nonzero block number
blkLen = 6;       % block length

SNR = 15;         % Signal-to-noise ratio
iterNum = 100;    % number of experiments


for it = 1: iterNum
    fprintf('\n\nRunning %d:\n',it);

    % Generate the known matrix with columns draw uniformly from the surface of a unit hypersphere
    Phi = randn(M,N);
    Phi = Phi./(ones(M,1)*sqrt(sum(Phi.^2)));
    
    % generate nonzero block coefficients
    
    beta = ones(blkNum,1)*0.99;
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

    IND = find(abs(Wgen)>0);
    
    
   % ======================================================================
   %             Algorithm Comparison
   % ======================================================================
   
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Benchmark
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    supt = find( abs(Wgen)>0);
    x_ls = pinv(Phi(:,supt)) * Y;
    mse_bench(it) = (norm(Wgen(supt) - x_ls,'fro')/norm(Wgen,'fro'))^2; 
    fprintf('   bench:  MSE: %g;\n', mean(mse_bench));
        
    
    
    %================== BSBL-EM exploiting intra-block correlation ==============
    blkStartLoc = [1:blkLen:N];
    learnlambda = 1;         % when SNR < 25dB, set learnlambda = 1;
                             % when SNR >=25dB, set learnlambda = 2;
                             % when no noise present, set learnlambda = 0;
    tic;
    
    Result1 = BSBL_EM(Phi,Y,blkStartLoc,learnlambda);

    t_em(it) = toc;
    mse_em(it) = (norm(Wgen - Result1.x,'fro')/norm(Wgen,'fro'))^2;


    fprintf('BSBL-EM(learn correlation) : time: %4.3f, MSE: %g\n',mean(t_em),mean(mse_em));
 
    
    
    
    %================== BSBL-EM ignoring intra-block correlation ==============

    tic;
    
    Result2 = BSBL_EM(Phi,Y,blkStartLoc,learnlambda,'learntype',0);

    t_em2(it) = toc;
    mse_em2(it) = (norm(Wgen - Result2.x,'fro')/norm(Wgen,'fro'))^2;


    fprintf('BSBL-EM(ignore correlation): time: %4.3f, MSE: %g\n\n',mean(t_em2),mean(mse_em2));
    
    
    
    
    %================== BSBL-BO exploiting intra-block correlation ==============
    
    learnlambda3 = 2;        % when SNR < 10dB, set learnlambda = 1;
                             % when SNR >=10dB, set learnlambda = 2;
                             % when no noise present, set learnlambda = 0;
    tic;
    
    Result3 = BSBL_BO(Phi,Y,blkStartLoc,learnlambda3);

    t_bo(it) = toc;
    mse_bo(it) = (norm(Wgen - Result3.x,'fro')/norm(Wgen,'fro'))^2;


    fprintf('BSBL-BO(learn correlation) : time: %4.3f, MSE: %g\n',mean(t_bo),mean(mse_bo));
    
    
    %================== BSBL-BO ignoring intra-block correlation ==============
    
    learnlambda3 = 2;        % when SNR < 10dB, set learnlambda = 1;
                             % when SNR >=10dB, set learnlambda = 2;
                             % when no noise present, set learnlambda = 0;
    tic;
    
    Result4 = BSBL_BO(Phi,Y,blkStartLoc,learnlambda3,'learnType',0);

    t_bo2(it) = toc;
    mse_bo2(it) = (norm(Wgen - Result4.x,'fro')/norm(Wgen,'fro'))^2;


    fprintf('BSBL-BO(ignore correlation) : time: %4.3f, MSE: %g\n',mean(t_bo2),mean(mse_bo2));
    
end