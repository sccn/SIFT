% In this demo we run 4 algorithms in a noisy experiment when the block partition is unknown.

clear; 


n = 500;     % signal length
N = 200;     % sensor measurement number
 
K = 50;      % sparsity (nonzero elements in the signal)
q = 6;       % group number (number of nonzero blocks)

SNR = 15;    % Signal-to-noise ratio
trial = 10;  % trial number

runBSBL_EM    = 1;    % run BSBL-EM (h=4 and h=8); you can turn off the algorithm by setting runBSBL_EM = 0;
runBSBL_BO    = 1;    % run BSBL-BO (h=4 and h=8)
runEBSBL_BO   = 1;    % run EBSBL-BO (h=4 and h=8)


for it = 1 : trial
    fprintf('\nTrial No.%d \n',it);
    
    % generate sparse signal
    x=zeros(n,1);
    r=abs(randn(q,1)); r=r+1; r=round(r*K/sum(r)); r(q)=K-sum(r(1:q-1));
    g=round(r*n/K); g(q)=n-sum(g(1:q-1));
    g_cum=cumsum(g);
    
    for i=1:q,
        
        % intra-block correlation
        beta = 1 - 0.05*rand;  
        
        % generate i-th group with intra-block correlation
        seg = [];
        seg(1) = randn;
        for kk = 2 : r(i)
            seg(kk,1) = beta*seg(kk-1,1) + sqrt(1-beta^2) * randn;
        end
        loc=randperm(g(i)-r(i)-2);
        x_tmp=zeros(g(i), 1);
        x_tmp(loc(1)+2:loc(1)+1+r(i))= seg; % generate coherent blocks
        
        x(g_cum(i)-g(i)+1:g_cum(i), 1)=x_tmp;
    end

    % true support
    supt = find( abs(x)>0);
    

    % generate the known matrix (random Gaussian matrix)
    Phi = randn(N,n);
    Phi = Phi./(ones(N,1)*sqrt(sum(Phi.^2)));

    % noiseless signal
    measure = Phi * x;

    % Observation noise
    stdnoise = std(measure)*10^(-SNR/20);
    noise = randn(N,1) * stdnoise;

    % Noisy signal
    Y = measure + noise;

    
    %============ benchmark (least square with given support) =============
     
    x_ls = pinv(Phi(:,supt)) * Y;
    mse_bench(it) = (norm(x(supt) - x_ls,'fro')/norm(x,'fro'))^2; 
    fprintf('benchmark    :              MSE: %g;\n', mean(mse_bench));
    %======================================================================
    

    
    %=========================================================
    %                BSBL-EM with h = 4
    %=========================================================
    if runBSBL_EM    == 1  

    h_em1 = 4;
    blkStartLoc1 = [1:h_em1:n];
    learnlambda_em1 = 1;          % when SNR is low (e.g. < 20dB), set to 1
                                  % when SNR is high (e.g. <20dB), set to 2
                                  % in noiseless cases, set to 0
    
    tic;
    Result1 = BSBL_EM(Phi, Y, blkStartLoc1, learnlambda_em1);
    
    t_em1(it) = toc;
    mse_em1(it) = (norm(x - Result1.x,'fro')/norm(x,'fro'))^2;


    fprintf('BSBL-EM (h=4): time: %4.3f, MSE: %g\n',mean(t_em1),mean(mse_em1));
    
    
    
    %=========================================================
    %                BSBL-EM with h = 8
    %=========================================================
    h_em2 = 8;
    blkStartLoc2 = [1:h_em2:n];
    learnlambda_em2 = 1;          % when SNR is low (e.g. < 20dB), set to 1
                                  % when SNR is high (e.g. <20dB), set to 2
                                  % in noiseless cases, set to 0
    tic;
    Result2 = BSBL_EM(Phi, Y, blkStartLoc2, learnlambda_em2);
    
    t_em2(it) = toc;
    mse_em2(it) = (norm(x - Result2.x,'fro')/norm(x,'fro'))^2;


    fprintf('BSBL-EM (h=8): time: %4.3f, MSE: %g\n',mean(t_em2),mean(mse_em2));
 
    end
    
 
    
    
    %=========================================================
    %                BSBL-BO with h = 4
    %=========================================================
    if runBSBL_BO    == 1  

    h_bo3 = 4;
    blkStartLoc3 = [1:h_bo3:n];
    learnlambda_bo3 = 2;    % set it to 2 if SNR > 10dB
    
    tic;
    Result3 = BSBL_BO(Phi, Y, blkStartLoc3, learnlambda_bo3);
    
    t_bo3(it) = toc;
    
    mse_bo3(it) = (norm(x - Result3.x,'fro')/norm(x,'fro'))^2;

 
    fprintf('BSBL-BO (h=4): time: %4.3f, MSE: %g\n',mean(t_bo3),mean(mse_bo3));
    
    
    
    %=========================================================
    %                BSBL-BO with h = 8
    %=========================================================
    h_bo4 = 8;
    blkStartLoc4 = [1:h_bo4:n];
    learnlambda_bo4 = 2;   % set it to 2 if SNR > 10dB
    
    tic;
    Result4 = BSBL_BO(Phi, Y, blkStartLoc4, learnlambda_bo4);
    
    t_bo4(it) = toc;
    
    mse_bo4(it) = (norm(x - Result4.x,'fro')/norm(x,'fro'))^2;

 
    fprintf('BSBL-BO (h=8): time: %4.3f, MSE: %g\n',mean(t_bo4),mean(mse_bo4));
    
    end
    
    
    %=========================================================
    %               EBSBL-BO with h = 4
    %=========================================================
    if runEBSBL_BO   == 1

    Len5 = 4;
    noiseFlag5 = 1;   % = 1: strong noise;
                      % = 2: small noise;  
                      % = 0 : no noise
    tic;
    Result5 = EBSBL_BO(Phi, Y, Len5, noiseFlag5);
    t_ecv5(it) = toc;
    
    mse_ecv5(it) = (norm(x - Result5.x,'fro')/norm(x,'fro'))^2;


    fprintf('EBSBL-BO(h=4): time: %4.3f, MSE: %g\n',mean(t_ecv5),mean(mse_ecv5));
    
    
 
    %=========================================================
    %               EBSBL-BO with h = 8
    %=========================================================
    Len6 = 8;
    noiseFlag6 = 1;   % = 1: strong noise;
                      % = 2: small noise;  
                      % = 0 : no noise
    tic;
    Result6 = EBSBL_BO(Phi, Y, Len6, noiseFlag6);
    t_ecv6(it) = toc;
    
    mse_ecv6(it) = (norm(x - Result6.x,'fro')/norm(x,'fro'))^2;


    fprintf('EBSBL-BO(h=8): time: %4.3f, MSE: %g\n\n',mean(t_ecv6),mean(mse_ecv6));
    
    end;
    


    
end


