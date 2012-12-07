function Result = BSBL_L1_noise(Phi, y, blkLen, iterNum)
% BSBL-L1: Recover block sparse signal (1D) exploiting intra-block correlation, given the block partition.
%          It uses the group Basis Pursuit to implement each iteration
%          The functions to perform group Basis Pursuit are downloaded from
%          the SGPL1 software: http://www.cs.ubc.ca/labs/scl/spgl1/
%
%    Note: This code is only suitable for noisy cases when each block
%          has the identitcal known block-size
%
% ============================== INPUTS ============================== 
%   Phi         : M X N known matrix
%
%   y           : M X 1 measurement vector 
%
%   blkLen      : length of each block (assuming identical block size)
%   
%   iterNum     : iteration number (the first iteration is just a standard
%                 group Basis Pursuit). 
%               [ Suggest setting: iterNum = 3 in most noisy cases; 
%                                  iterNum = 6 in strongly noisy cases (SNR<5dB)]
%                 
%
% ==============================  OUTPUTS ============================== 
%   Result : A structured data with:
%     Result.x     : a solution matrix; the i-th column is the solution at
%                    the i-th iteration
%     Result.time  : a vector; the i-th entry records the time at the i-th
%                    iteration (note: not the used time at the i-th iteration)
%     Result.B     : the final B matrix
%
%
% ================================= See Also =============================
%   BSBL_BO,   BSBL_EM,    TMSBL,    TSBL      
%
% ================================ Reference =============================
%   [1] Zhilin Zhang, Bhaskar D. Rao, Extension of SBL Algorithms for the 
%       Recovery of Block Sparse Signals with Intra-Block Correlation, 
%       available at: http://arxiv.org/abs/1201.0862
%
%   [2] For more info:  http://dsp.ucsd.edu/~zhilin/BSBL.html
%
% ============= Author =============
%   Zhilin Zhang (z4zhang@ucsd.edu)
%
% ============= Version =============
%   1.1 (01/24/2012)
%   1.0 (08/27/2011)
%


[M,N] = size(Phi);            % size 

g = N/blkLen;                 % block number


% Initialize B
% r = 0.0;
% for ip = 1 : blkLen
%     bs(ip) = r^(ip-1);
% end
B = eye(blkLen);
invB = B;
 

X_hat = zeros(N,iterNum);       % used to record the solution at each iteration

PRUNE_GAMMA = 1e-3;             % threshold to prune out small gamma_i

lambda0 = (std(y)*1e-2)^2;      % this parameter is used for group Basis Pursuit

time0 = cputime;                % record the starting time

% Iterative reweighting
for iter = 1 : iterNum
    
    if iter == 1               
        wgn = ones(g,1);
        gamma = ones(g,1); 
        keep_list = [1:g]';
        usedNum = length(keep_list);
        x = zeros(N,1);
        
    else
        
        % =============== computing new gamma ===============
        for iq = 1 : usedNum
            seg = (iq-1)*blkLen+1 : iq*blkLen;
            gamma(iq) = sqrt(x(seg)'* invB * x(seg))/wgn(iq);
            
            if gamma(iq) < PRUNE_GAMMA
                gamma(iq) = 0;
            end;
        end
          
        
        % =============== computing lambda ==================
        lambda = norm(y - Phi*x)^2/M;   % this parameter is used for the updating of B, gamma, and weighting
        
        
%         % ================== computing B ====================
%         B = zeros(blkLen);
%         for gidx = 1 : usedNum
%             blkInd = (gidx-1)*blkLen + 1 : gidx*blkLen;
%             
%             if gamma(gidx) >= PRUNE_GAMMA
%             B = B + x(blkInd) * x(blkInd)'/gamma(gidx);
%             end
%         end
%         b0 = mean(diag(B)); b1 = mean(diag(B,1));
%         if abs(b1/b0) > 0.99, 
%             r = sign(b1/b0) * 0.99;
%         else
%             r = b1/b0;
%         end
%         for db = 1 : blkLen
%             ars(db) = r^(db-1);
%         end
%          
%         B = toeplitz(ars);   
% 
%         invB = inv(B);
        
        
        % =============== computing new weights ================
        PGP = lambda * eye(M) + Phi*kron(diag(gamma),B)*Phi';
        invPGP = pinv(PGP);
        wgn = [];
        for iq = 1 : usedNum
            seg = (iq-1)*blkLen+1 : iq*blkLen;
            wgn(iq) = sqrt(trace(B*Phi(:,seg)'*invPGP*Phi(:,seg)));
            
        end
    end

    
    % transform to a standard group Lasso type problem
    invR = kron( diag(1./wgn.^2),sqrtm(B));
    Phi_tild = Phi * invR;
    
    
    % ==========================================================================
    % Use the standard group Basis Pursuit to solve the group Lasso type problem 
    %===========================================================================
    x_old = x;
    
    % -- calcualte the 1st parameter for group Basis Pursuit
    Gsize = blkLen; nGroups = N/Gsize; groups = [1:nGroups];
    groups = repmat(groups, [Gsize, 1]);
    groups = groups(:);
     
    % -- calcualte the 2nd parameter for group Basis Pursuit
    sigma = sqrt(lambda0) * sqrt(M+2*sqrt(2*M));
    
    % -- calcualte the 3rd parameter for group Basis Pursuit
    opts = spgSetParms('verbosity',0,'iterations',100,'bpTol',1e-04,'optTol',1e-04);
      
    % -- run group basis pursuit
    x = spg_group(Phi_tild,y,groups,sigma,opts);
     
    % transform the solution to the original model
    x = invR * x;
    
    % record the past time at this iteration
    elapseTime(iter,1) = cputime - time0;
    
    
    %===========================================================================
    %     post processing 
    %===========================================================================
    x0 = zeros(N,1); activeSet = [];
    for kk = 1 : usedNum
        activeSet = [activeSet, (keep_list(kk)-1)*blkLen+1:keep_list(kk)*blkLen ];
    end
    x0(activeSet) = x;
    X_hat(:,iter) = x0;

    
    %===========================================================================
    %    decide if the algorithm has converged
    %===========================================================================
    if norm(x-x_old)/norm(x) < 1e-2    
        
        for jj = iter + 1 : iterNum
            X_hat(:,jj) = x0;
            elapseTime(jj,1) = cputime - time0;
        end
        break;
    end
end

% Output
Result.x = X_hat;
Result.time = elapseTime;
Result.B = B;





