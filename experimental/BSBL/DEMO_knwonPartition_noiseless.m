% This demo shows how to set input parameters of BSBL-EM and BSBL-BO 
% for noiseless experiments when block partition is known.
%
% Also it shows the performance improvement if exploiting intra-block correlation

clear; 

cor = 0.9;         % true intra-block correlation

% problem dimension
M = 100;          % row number of the matrix 
N = 300;          % column number

blkNum = 20;       % number of nonzero blocks
blkLen = 4;        % block length

iterNum = 20;       % number of experiments
 

        
for it = 1: iterNum
    fprintf('\n\nRunning %d:\n',it);

    % Generate the known matrix with columns draw uniformly from the surface of a unit hypersphere
    Phi = randn(M,N);
    Phi = Phi./(ones(M,1)*sqrt(sum(Phi.^2)));
 
    % generate nonzero block coefficients
    beta = ones(blkNum,1) * cor; 
    blks(:,1) = randn(blkNum,1);
    for i = 2 : blkLen
        blks(:,i) = beta .* blks(:,i-1) + sqrt(1-beta.^2).*(ones(blkNum,1).*randn(blkNum,1));
    end
    blks = blks./( sqrt(sum(blks.^2,2)) * ones(1,blkLen) );
    
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
    Y = Phi * Wgen;
    

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             BSBL-EM exploiting intra-block correlation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    groupStartLoc = 1:blkLen:N;         % start location of each block
    
    % Learn B
    learnlambda = 0;
    tic;
    Result1 = BSBL_EM(Phi,Y,groupStartLoc, learnlambda);
    t_csbl1(it) = toc;
    
    % mse
    mse_cSBL1(it) = (norm(Wgen - Result1.x,'fro')/norm(Wgen,'fro'))^2;
    
    % calculate failure rate
    if mse_cSBL1(it) > 1e-6
        fail_cSBL1(it) = 1;
    else
        fail_cSBL1(it) = 0;
    end
    
    
    fprintf('BSBL-EM(learn B ): time: %4.3f, MSE: %g, Success: %4.3f,count: %d; \n',...
        mean(t_csbl1),mean(mse_cSBL1),1-mean(fail_cSBL1),Result1.count);
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             BSBL-EM when ignores correlation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    learnlambda2 = 0;
    tic;
    Result2 = BSBL_EM(Phi,Y,groupStartLoc,learnlambda2, 'LEARNTYPE', 0);   
    t_csbl2(it) = toc;
    
    % mse
    mse_cSBL2(it) = (norm(Wgen - Result2.x,'fro')/norm(Wgen,'fro'))^2;
    
    % calculate failure rate
    if mse_cSBL2(it) > 1e-6
        fail_cSBL2(it) = 1;
    else
        fail_cSBL2(it) = 0;
    end
    
    fprintf('BSBL-EM(ignore B): time: %4.3f, MSE: %g, Success: %4.3f,count: %d; \n',...
        mean(t_csbl2),mean(mse_cSBL2),1-mean(fail_cSBL2),Result2.count);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              BSBL-BO exploiting intra-block correlation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Learn B
    learnlambda3 = 0;
    tic;
    Result3 = BSBL_BO(Phi,Y,groupStartLoc, learnlambda3);
    t_csbl3(it) = toc;
    
    % mse
    mse_cSBL3(it) = (norm(Wgen - Result3.x,'fro')/norm(Wgen,'fro'))^2;
    
    % calculate failure rate
    if mse_cSBL3(it) > 1e-6
        fail_cSBL3(it) = 1;
    else
        fail_cSBL3(it) = 0;
    end
    
    
    fprintf('BSBL-BO(learn B ): time: %4.3f, MSE: %g, Success: %4.3f,count: %d; \n',...
        mean(t_csbl3),mean(mse_cSBL3),1-mean(fail_cSBL3),Result3.count);
    
    
 

end


