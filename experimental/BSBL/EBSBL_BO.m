function Result = EBSBL_BO(Phi, y, Len, noiseFlag, varargin)

% EBSBL-BO: Recover block sparse signal (1D) exploiting intra-block correlation, with unknown block partition.
% 
% ============================== INPUTS ============================== 
%   Phi         : N X M known matrix
%
%   y           : N X 1 measurement vector 
%
%   Len         : the block size in the expanded model, i.e. the parameter h in the paper [1]
%   
%   noiseFlag   :  noiseFlag == 0, noiseless experiment
%                  noiseFlag == 1, strongly noisy experiment (eg. SNR < 20dB)
%                  noiseFlag == 2, lightly noisy experiment (eg. SNR >20dB)
%
% [varargin values -- in most cases you can use the default values]
%
%  'LEARNTYPE'   : LEARNTYPE = 0: Ignore intra-block correlation
%                  LEARNTYPE = 1: Exploit intra-block correlation 
%                 [ Default: LEARNTYPE = 1 ]
%
%  'Lambda'      : lambda value
%                 [ Default: 10^{-14} in noiseless; std(y)*1e-2 in noisy case]
%
%  'PRUNE_GAMMA'  : threshold to prune out small gamma_i 
%                   (generally, 10^{-3} or 10^{-2})
%
%  'MAX_ITERS'    : Maximum number of iterations.
%                 [ Default value: MAX_ITERS = 500 ]
%
%  'EPSILON'      : Solution accurancy tolerance parameter 
%                 [ Default value: EPSILON = 1e-7   ]
%
%  'PRINT'        : Display flag. If = 1: show output; If = 0: supress output
%                 [ Default value: PRINT = 0        ]
%
% ==============================  OUTPUTS ============================== 
%   Result : A structured data with:
%      Result.x          : the estimated block sparse signal
%      Result.blkind     : indexes of nonzero blocks in the expanded sparse block model
%      Result.gamma_est  : the gamma values associated with all the blocks in the model
%      Result.B          : the final value of the matrix B
%      Result.count      : iteration times
%      Result.lambda     : the final value of lambda
%
%
% ========================= Command examples  =============================
%   < Often-used command >
%    For noisy environment:
%          
%           Result =  EBSBL_BO(Phi, y, Len, 1);  
%
%    For noiseless cases:
%          
%           Result =  EBSBL_BO(Phi, y, Len, 0);  
%
%    < Full command >
%           Result =  EBSBL_BO(Phi, y, Len, ...
%                                      'LEARNTYPE', 1,...
%                                      'Lambda', 1e-3,...
%                                      'MAX_ITERS', 800,...
%                                      'EPSILON', 1e-8,...
%                                      'PRINT',0);
%
% ================================= See Also =============================
%   BSBL_BO,   BSBL_EM,  BSBL_L1,  EBSBL_L1,  TMSBL,    TSBL      
%
% ================================ Reference =============================
%   [1] Zhilin Zhang, Bhaskar D. Rao, Extension of SBL Algorithms for the 
%       Recovery of Block Sparse Signals with Intra-Block Correlation, 
%       available at: http://arxiv.org/abs/1201.0862
%
%   [2] For more info:  http://dsp.ucsd.edu/~zhilin/BSBL.html
%
%   [3] Research Blog:  http://marchonscience.blogspot.com/
%
% ============= Author =============
%   Zhilin Zhang (z4zhang@ucsd.edu)
%
% ============= Version =============
%   1.1 (05/28/2012)
%   1.0 (01/22/2012)
%


% scaling...
scl = std(y);
if (scl < 0.4) | (scl > 1)
    y = y/scl*0.4;
end


% Dimension of the Original Problem
[N M] = size(Phi); 

% generalize the block structure
C0 = sqrtm(eye(Len));
p = M - Len + 1;
A = zeros(N,p*Len);
for i = 1 : p
    Ci = zeros(M,Len);
    Ci( [i:i+Len-1], : ) = C0;
    C{i} = Ci;
    A(:,(i-1)*Len+1:i*Len) = Phi * C{i};
end


% Default Control Parameters ***
EPSILON     = 1e-7;       % solution accurancy tolerance
MAX_ITERS   = 500;        % maximum iterations
PRINT       = 0;          % don't show progress information
LEARNTYPE     = 1;        % adaptively estimate the covariance matrix B

if noiseFlag == 0  
    lambda = 1e-12;   
    PRUNE_GAMMA = 1e-3;
elseif noiseFlag == 1  
    lambda = std(y) * 1e-2;    
    PRUNE_GAMMA = 1e-2;
elseif noiseFlag == 2
    lambda = std(y) * 1e-3;
    PRUNE_GAMMA = 1e-2;
else
    error(['Unrecognized Value for Input Argument ''noiseFlag''']);
end

 

if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'learntype'
                LEARNTYPE = varargin{i+1};
                if LEARNTYPE ~= 1 & LEARNTYPE ~= 0
                    error(['Unrecognized Value for Input Argument ''LEARNTYPE''']);
                end
            case 'lambda'
                lambda = varargin{i+1};
            case 'prune_gamma'
                PRUNE_GAMMA = varargin{i+1}; 
            case 'epsilon'   
                EPSILON = varargin{i+1}; 
            case 'print'    
                PRINT = varargin{i+1}; 
            case 'max_iters'
                MAX_ITERS = varargin{i+1};  
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end


if PRINT
    fprintf('\n====================================================\n');
    fprintf('           Running EBSBL-BO ....... \n');
    fprintf('           Information about parameters...\n');
    fprintf('====================================================\n');
    fprintf('PRUNE_GAMMA  : %e\n',PRUNE_GAMMA);
    fprintf('lambda       : %e\n',lambda);  
    fprintf('LearnType    : %d\n',LEARNTYPE);
    fprintf('EPSILON      : %e\n',EPSILON);
    fprintf('MAX_ITERS    : %d\n\n',MAX_ITERS);
end



%% Initialization 
if LEARNTYPE == 1,
    for i = 1 : Len, b(i) = 0.9^(i-1); end
    B = toeplitz(b);
else
    B = eye(Len);
end

gamma = ones(p,1);
keep_list = [1:p]';
usedNum = length(keep_list);
mu_x = zeros(p*Len,1);
count = 0;


%% Iteration
while (1)
    count = count + 1;

    %=========== Prune weights as their hyperparameters go to zero ==============
    if (min(gamma) < PRUNE_GAMMA)
        index = find(gamma > PRUNE_GAMMA);
        usedNum = length(index);
        
        % prune gamma 
        gamma = gamma(index);  
         
        % Find corresponding columns 
        extIndex = zeros(usedNum,Len);
        for i = 1 : usedNum
            extIndex(i,:) = (index(i)-1)*Len+1 : index(i)*Len;
        end
        extIndex = reshape(extIndex',usedNum*Len,1);
        A = A(:,extIndex);
          
        % record the remained indexes
        keep_list = keep_list(index);
    end
    
    %=================== Compute new weights =================
    mu_old = mu_x;
    
    DBD = zeros(N);
    for i = 1 : usedNum
        DBD = DBD + A(:, (i-1)*Len+1: i*Len ) * gamma(i) * B * A(:, (i-1)*Len+1: i*Len )';
    end
    H = A'/(DBD + lambda * eye(N));
    Ht = H*y;      HD = H * A;
    
    mu_x = zeros(usedNum*Len,1);
    Sigma_x = repmat(zeros(Len),[1 1 usedNum]);
    gamma_old = gamma;  
    Cov_x = Sigma_x;
    Bnew = zeros(Len);
    for i = 1 : usedNum
        seg = [(i-1)*Len+1 : i*Len];
        mu_x(seg) = gamma(i) * B * Ht(seg);       % solution
        Sigma_x(:,:,i) = gamma(i) * B - gamma(i)^2 * B * HD(seg,seg) * B; 
        Cov_x(:,:,i) = Sigma_x(:,:,i) + mu_x(seg) * mu_x(seg)';
        
        % Update gamma
        gamma(i) = gamma_old(i)*norm( sqrtm(B)*Ht(seg) )/sqrt(trace(HD(seg,seg)*B));
        
        % Estimate covariance matrix
        if LEARNTYPE == 1,
            Bnew = Bnew + Cov_x(:,:,i)/gamma(i); 
        end
        
    end
    
    % reconstruct B_i for each block
    if LEARNTYPE == 1,
        b = (mean(diag(Bnew,1))/mean(diag(Bnew)));
        if abs(b) > 0.99, b = 0.99*sign(b); end;
        bs = [];
        for j = 1 : Len, bs(j) = (b)^(j-1); end;
        B = toeplitz(bs);
    else
        B = eye(Len);
    end
    
%     % Estimate lambda
%     if LearnLambda == 1
%         lambdaComp = 0; 
%         for i = 1 : usedNum
%             currentSeg = [(i-1)*Len+1 : i*Len];
%             lambdaComp = lambdaComp + trace(A(:,currentSeg)*Sigma_x(:,:,i)*A(:,currentSeg)');
%         end
%         lambda = norm(y - A * mu_x,2)^2/N + lambdaComp/N; 
%         
%     elseif LearnLambda == 2
%         lambdaComp = 0;
%         invB = inv(B);
%         for i = 1 : usedNum
%             lambdaComp = lambdaComp + trace(Sigma_x(:,:,i)*invB)/gamma_old(i);
%         end
%         lambda = norm(y - A * mu_x,2)^2/N + lambdaComp/N; 
%     end
    

    

    % ================= Check stopping conditions, etc. ==============
    if (PRINT) disp([' iters: ',num2str(count),'   num coeffs: ',num2str(usedNum), ...
            '   gamma change: ',num2str(max(abs(gamma - gamma_old)))]); end;
    if (count >= MAX_ITERS), if PRINT, fprintf('Reach max iterations. Stop\n\n'); end; break;  end;

    if (size(mu_x) == size(mu_old))
        dmu = max(max(abs(mu_old - mu_x)));
        if (dmu < EPSILON)  break;  end;
    end;
    

end;



%% Reconstruct the signal

x_est = zeros(M,1);
for i = 1 : usedNum
    x_est = x_est + C{keep_list(i)} * mu_x( (i-1)*Len+1 : i*Len );
end

%% output results

if (scl < 0.4) | (scl > 1)
    Result.x = x_est * scl/0.4;
else
    Result.x = x_est;
end

Result.blkind = sort(keep_list);
gamma_est = zeros(p,1); gamma_est(keep_list,1) = gamma;  
Result.gamma_est = gamma_est;
Result.B = B;
Result.count = count;
Result.lambda = lambda;
return;
