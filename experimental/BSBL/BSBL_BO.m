function Result = BSBL_BO(varargin)

% BSBL-BO: Recover block sparse signal (1D) exploiting intra-block correlation, given the block partition.
%
%          The algorithm solves the inverse problem for the block sparse
%          model with known block partition:
%                        y = Phi * x + v
%
%
% ============================== INPUTS ============================== 
%   Phi         : N X M known matrix
%
%   y           : N X 1 measurement vector 
%
%   blkStartLoc : Start location of each block
%   
%   LearnLambda : (1) If LearnLambda = 1, use the lambda learning rule for very LOW SNR cases (SNR<10dB)
%                     (using lambda=std(y)*1e-2 or user-input value as initial value)
%                 (2) If LearnLambda = 2, use the lambda learning rule for medium noisy cases (SNR>10dB) 
%                     (using lambda=std(y)*1e-2 or user-input value as initial value)
%                 (3) If LearnLambda = 0, do not use the lambda learning rule 
%                     ((using lambda=1e-14 or user-input value as initial value)
%                 
%
% [varargin values -- in most cases you can use the default values]
%
%   'LEARNTYPE'  : LEARNTYPE = 0: Ignore intra-block correlation
%                  LEARNTYPE = 1: Exploit intra-block correlation 
%                 [ Default: LEARNTYPE = 1 ]
%
%  'PRUNE_GAMMA'  : threshold to prune out small gamma_i 
%                   (generally, 10^{-3} or 10^{-2})
%
%  'LAMBDA'       : user-input value for lambda
%                  [ Default: LAMBDA=1e-14 when LearnLambda=0; LAMBDA=std(y)*1e-2 in noisy cases]
%
%  'MAX_ITERS'    : Maximum number of iterations.
%                 [ Default value: MAX_ITERS = 600 ]
%
%  'EPSILON'      : Solution accurancy tolerance parameter 
%                 [ Default value: EPSILON = 1e-8   ]
%
%  'PRINT'        : Display flag. If = 1: show output; If = 0: supress output
%                 [ Default value: PRINT = 0        ]
%
% ==============================  OUTPUTS ============================== 
%   Result : 
%      Result.x          : the estimated block sparse signal
%      Result.gamma_used : indexes of nonzero groups in the sparse signal
%      Result.gamma_est  : the gamma values of all the groups of the signal
%      Result.B          : the final value of the B
%      Result.count      : iteration times
%      Result.lambda     : the final value of lambda
%
%
% ========================= Command examples  =============================
%   < Often-used command >
%    For most noisy environment (SNR > 10dB):
%          
%           Result =  BSBL_BO(Phi, y, blkStartLoc, 2);  
%
%    For very low SNR cases (SNR < 10 dB):
%           
%           Result =  BSBL_BO(Phi, y, blkStartLoc, 1);   
%
%    For noiseless cases:
%          
%           Result =  BSBL_BO(Phi, y, blkStartLoc, 0);  
%
%    To recover non-Sparse structured signals (noiseless):
%           Result = BSBL_BO(Phi,y,groupStartLoc,0,'prune_gamma',-1);
%                    ('prune_gamma' can be set any non positive constant)
%
%    < Full-Command Example >
%           Result =  BSBL_BO(Phi, y, blkStartLoc, learnlambda, ...
%                                                 'LEARNTYPE', 1,...
%                                                 'PRUNE_GAMMA',1e-2,...
%                                                 'LAMBDA',1e-3,...
%                                                 'MAX_ITERS', 800,...
%                                                 'EPSILON', 1e-8,...
%                                                 'PRINT',0);
%
% ================================= See Also =============================
%   EBSBL_BO,   BSBL_EM,  BSBL_L1,  EBSBL_L1,  TMSBL,    TSBL      
%
% ================================ Reference =============================
%   [1] Zhilin Zhang, Bhaskar D. Rao, Extension of SBL Algorithms for the 
%       Recovery of Block Sparse Signals with Intra-Block Correlation, 
%       available at: http://arxiv.org/abs/1201.0862
%  
%   [2] Zhilin Zhang, Tzyy-Ping Jung, Scott Makeig, Bhaskar D. Rao, 
%       Low Energy Wireless Body-Area Networks for Fetal ECG Telemonitoring 
%       via the Framework of Block Sparse Bayesian Learning, 
%       available at: http://arxiv.org/pdf/1205.1287v1.pdf
%
%   [3] webpage: http://dsp.ucsd.edu/~zhilin/BSBL.html, or
%                https://sites.google.com/site/researchbyzhang/bsbl
%
% ============= Author =============
%   Zhilin Zhang (z4zhang@ucsd.edu, zhangzlacademy@gmail.com)
%
% ============= Version =============
%   1.4 (07/23/2012) debug
%   1.3 (05/30/2012) make faster
%   1.2 (05/28/2012)
%   1.1 (01/22/2012)
%   1.0 (08/27/2011)
%

arg_define(varargin, ...
    arg_norep({'Phi','DesignMatrix','A'},mandatory,[],'The design matrix. This is the data matrix (ie X).'), ...
    arg_norep({'y','TargetVector'},mandatory,[],'The target vector'), ...
    arg_norep({'blks','Blocks'},mandatory,[],'Block (Group) size. Array containing the number of elements in each block (group) of coefficients. ie: blks = [5 5 10] is two blocks of size 5 followed by a block of size 10.'), ...
    arg_nogui({'designMatrixBlockSize','DesignMatrixBlockSize'},[],[],'Design matrix structure. If empty, do not exploit structure for faster computation. If non-empty, the design matrix consists of identical blocks along the main diagonal. ''DesignMatrixBlockSize'' is the [numrow numcol] size of each block.'), ...
    arg({'LearnLambda','LambdaLearningRule'},'LowNoise',{'off','LowNoise','HighNoise'},{'Lambda learning rule',sprintf(['\n' ...
                                                                                          '''LowNoise'' :  Use the lambda learning rule for very LOW SNR cases (SNR<10dB). Uses lambda=std(y)*1e-2 or user-input value as initial value.\n\n' ...
                                                                                          '''HighNoise'':  Use the lambda learning rule for medium noisy cases (SNR>10dB). Uses lambda=std(y)*1e-2 or user-input value as initial value.\n\n' ...
                                                                                          '''off'':        Do not use the lambda learning rule. Uses lambda=1e-14 or user-input value as initial value' ...
                                                                                          ])}), ...
    arg({'PRUNE_GAMMA','GammaPruningThresh'},[],[0 Inf],{'Gamma Pruning threshold',sprintf(['\n' ...
                                                           'Threshold for prunning small hyperparameters gamma_i. Blocks of parameters are pruned when their ''power'' (gamma_i) is small.\n'    ...
                                                           'In noisy cases, you can set PRUNE_GAMMA = 1e-3 or 1e-4. \n'   ...
                                                           'In strong noisy cases (e.g. SNR <= 6 dB), set PRUNE_GAMMA = 1e-2 for better performance.' ...
                                                         ])}), ...
    arg({'LAMBDA','InitialLambda','Lambda'},[],[],'Initial regularization value (lambda). If empty, the initial value is estimated from the data.','type','denserealdouble'), ...                                                        
    arg({'EPSILON','StoppingTolerance'},1e-8,[0 Inf],'Stopping tolerance. The algorithm terminates when the L2 norm of the change in parameter estimates is small than this value.'), ...
    arg({'LEARNTYPE','LearnCorrelation'},true,[],'Exploit correlation structure.'), ...
    arg({'MAX_ITERS','MaxIterations'},600,[1 Inf],'Maximum iterations'), ...
    arg({'PRINT','VerboseOutput','verb'},false,[],'Verbosity'), ...
    arg_nogui({'initState','InitialState'},[],[],'Initial BSBL state object for warm start.') ...
    );

% Error checks
if isstruct(initState) && ~isempty(initState) && ~all(isfield(initState,{'x' 'B'}))
	error('Invalid format for Input Argument ''initState''. This should be a structure such as output by BSBL_BO.');
end

% Modify some arguments
switch LearnLambda
    case 'off',         LearnLambda = 0;
    case 'LowNoise',    LearnLambda = 1;
    case 'HighNoise',   LearnLambda = 2;
end

% Determine scaling factor
scl = std(y);
if (scl < 0.4) || (scl > 1)
    y = y/scl*0.4;
end

% Specify some default values
if isempty(LAMBDA)
    LAMBDA = fastif(LearnLambda,scl*1e-2,1e-12);
end
if isempty(PRUNE_GAMMA)
    PRUNE_GAMMA = fastif(LearnLambda,1e-2,1e-3);
end

% Print some information
if PRINT
    fprintf('\n====================================================\n');
    fprintf('           Running BSBL-BO                \n');
    fprintf('           Information about parameters   \n');
    fprintf('====================================================\n');
    fprintf('PRUNE_GAMMA  : %e\n',PRUNE_GAMMA);
    fprintf('LAMBDA       : %e\n',LAMBDA);
    fprintf('LearnLambda  : %d\n',LearnLambda);    
    fprintf('LearnType    : %d\n',LEARNTYPE);
    fprintf('EPSILON      : %e\n',EPSILON);
    fprintf('MAX_ITERS    : %d\n\n',MAX_ITERS);
end


%% Initialization
[N,M] = size(Phi);
blkStartLoc = [1 cumsum(blks(1:length(blks)-1))+1];
Phi0 = Phi;
blkStartLoc0 = blkStartLoc;
blks0 = blks;
nblks = length(blkStartLoc);   % number of blocks
maxLen = max(blks);
if all(blks==blks(1))
    equalSize = 1;
else
    equalSize = 0;
end

% initialize Sigma0 (cell array containing each block's covariance matrix)
Sigma0 = arrayfun(@(x)speye(x), blks, 'UniformOutput',false);

if isempty(initState)
    % initialize to default state
    gamma       = ones(nblks,1);
    keep_list   = (1:nblks)';
    usedNum     = length(keep_list);
    mu_x        = zeros(M,1);
    count       = 0;
else
    % initialize based on input state
    gamma       = initState.gamma_est;
    keep_list   = (1:length(gamma))'; %initState.gamma_used;
    usedNum     = length(keep_list);
    mu_x        = initState.x;
    count       = 0;
    LAMBDA      = initState.lambda;
end

%      Result.x          : the estimated block sparse signal
%      Result.gamma_used : indexes of nonzero groups in the sparse signal
%      Result.gamma_est  : the gamma values of all the groups of the signal
%      Result.B          : the final value of the B
%      Result.count      : iteration times
%      Result.lambda     : the final value of lambda


%% Main Iteration Loop
while count < MAX_ITERS

    
    %===== Prune weights and blocks as their hyperparameters go to zero ===
    if (min(gamma) < PRUNE_GAMMA)
        blks2keep = find(gamma > PRUNE_GAMMA);
        usedNum   = length(blks2keep);
        keep_list = keep_list(blks2keep); 
        if isempty(keep_list), 
            fprintf('\n====================================================================================\n');
            fprintf('x becomes zero vector. The solution may be incorrect. \n');
            fprintf('Current ''prune_gamma'' = %g, and Current ''EPSILON'' = %g.\n',PRUNE_GAMMA,EPSILON);
            fprintf('Try smaller values of ''prune_gamma'' and ''EPSILON'' or normalize ''y'' to unit norm.\n');
            fprintf('====================================================================================\n\n');
            break; 
        end
        
        % Prune gamma and associated components in Sigma0
        % (retain only blocks with nonzero gamma)
        gamma  = gamma(blks2keep);
        Sigma0 = Sigma0(blks2keep);

        % Prune the design matrix Phi 
        % (retain only blocks with nonzero gamma)
        
        % get linear indices into blocks
        blkStartLoc = blkStartLoc0(blks2keep);
        blks        = blks0(blks2keep);
        blkInd2keep = zeros(1,sum(blks));
        curidx = 1;
        for k = 1:usedNum
            blkInd2keep(curidx:curidx+blks(k)-1) = blkStartLoc(k):blkStartLoc(k)+blks(k)-1;
            curidx = curidx+blks(k);
        end
        Phi = Phi0(:,blkInd2keep);
        
%         temp = [];
%         for k = 1 : usedNum
%             temp = [temp, Phi0(:,blkStartLoc(k):blkStartLoc(k)+blkLen{k}-1)];
%         end
%         Phi = temp;
%         clear temp;
    end

    % get linear indices of elements in each block
    blkLen = cellfun(@(x)size(x,1), Sigma0);
    blkInd = cell(1,usedNum);
    for gidx = 1:usedNum
        blkInd{gidx} = (gidx-1)*blkLen(gidx)+1 : gidx*blkLen(gidx);
    end
    
    %=================== Compute new weights =================
    mu_old = mu_x;
    
%     PhiBPhi = cellfun(...
%         @(Sigma0_i,blkInd_i) Phi(:,blkInd_i)*Sigma0_i*Phi(:,blkInd_i)', ...
%         Sigma0, blkInd,'UniformOutput',false);

    PhiBPhi = zeros(N);
    for i=1:usedNum
        PhiBPhi = PhiBPhi + Phi(:,blkInd{i})*Sigma0{i}*Phi(:,blkInd{i})';
    end
    
    H    = Phi'/(PhiBPhi + LAMBDA * speye(N));
    Hy   = H * y;      
    HPhi = H * Phi;
    
    % Compute mean and covariance of parameters
    mu_x    = cellfun(...
        @(Sigma0_i,blkInd_i) Sigma0_i*Hy(blkInd_i), ...
        Sigma0,blkInd,'UniformOutput',false);
    Sigma_x = cellfun( ...
        @(Sigma0_i,blkInd_i) Sigma0_i-Sigma0_i*HPhi(blkInd_i,blkInd_i)*Sigma0_i, ...
        Sigma0, blkInd,'UniformOutput',false);
    Cov_x   = cellfun( ...
        @(Sigma_xi,mu_xi) Sigma_xi + mu_xi*mu_xi', ...
        Sigma_x, mu_x,'UniformOutput',false);
    mu_x    = cell2mat(cellfun(@transpose,mu_x,'UniformOutput',false))';

    %=========== Learn correlation structure ===========
    if LEARNTYPE == 0
        % assume identity correlation structure
        [B invB] = deal(cellfun(@eye,blkLen,'UniformOutput',false));
       
    elseif LEARNTYPE == 1
        % constrain all the blocks to have the same correlation structure
        r0 = 0;
        r1 = 0;
        B0 = zeros(maxLen);
        
        for i = 1 : usedNum
            if equalSize == 0
                if blkLen(i) > 1
                    tmp = Cov_x{i}/gamma(i);
                    r0  = r0 + mean(diag(tmp));
                    r1  = r1 + mean(diag(tmp,1));
                end
            elseif equalSize == 1
                B0 = B0 + Cov_x{i}/gamma(i);
            end
        end

    end % if LEARNTYPE
        
    %%
    
    %=========== Learn correlation structure in blocks with Constraint 1 ===========
    % If blocks have the same size
    if (LEARNTYPE == 1) && (equalSize == 1)

        % Constrain all the blocks have the same correlation structure
        % (an effective strategy to avoid overfitting)
        b = (mean(diag(B0,1))/mean(diag(B0)));
        if abs(b) >= 0.99
            b = 0.99*sign(b); 
        end
        B0   = toeplitz(b.^((1:maxLen)-1));
        B    = repmat({B0},1,usedNum);
        invB = repmat({inverse(B0)},1,usedNum);
        
    
    % if blocks have different sizes
    elseif (LEARNTYPE == 1) && (equalSize == 0)
        r = r1/r0; 
        if abs(r) >= 0.99 
            r = 0.99*sign(r); 
        end
        B    = cellfun(@(Cov_x_i) toeplitz(r.^((1:size(Cov_x_i,1))-1)), ...
                       Cov_x, 'UniformOutput',false);
        invB = cellfun(@inverse,B,'UniformOutput',false);

    end

    % update lambda and each block's gamma
    % for speed, we use cellfun() throughout instead of looping over blocks
    % ---------------------------------------------------------------------
    gamma_old  = gamma;

    % get the indices into each block
    if LearnLambda
        blkLen = cellfun(@(x)size(x,1), Sigma_x);
    else
        blkLen = cellfun(@(x)size(x,1), Sigma0);
    end
    blkInd = cell(1,usedNum);
    for gidx = 1:usedNum
        blkInd{gidx} = (gidx-1)*blkLen(gidx)+1 : gidx*blkLen(gidx);
    end
    
    % compute gamma update factors...
    gamfun = @(B_i,blkInd_i) ...
               norm(sqrtm(B_i)*Hy(blkInd_i)) ...
               /sqrt(trace(HPhi(blkInd_i,blkInd_i)*B_i));
    gamupd = vec(cellfun(gamfun,B,blkInd));
    % ...update the gammas
    gamma  = gamma_old.*gamupd; 

    % compute weighted covariance Sigma0{i} = B{i}*gamma(i)
    Sigma0 = cellfun(@mtimes,B,arrayfun(@(x){x},gamma'),'UniformOutput',false);

    % update lambda using the appropriate learning rule
    switch LearnLambda
        case 1
            lambdafun  = @(Sigma_xi,blkInd_i) ...
                           trace( Phi(:,blkInd_i)*Sigma_xi*Phi(:,blkInd_i)' );
            lambdaComp = sum( cellfun(lambdafun,Sigma_x,blkInd) );
            LAMBDA     = norm(y - Phi * mu_x, 2)^2/N + lambdaComp/N; 
        case 2
            lambdafun  = @(Sigma_xi,invB_i,gamma_old_i) ...
                           trace( Sigma_xi*invB_i )/gamma_old_i;
            lambdatmp  = cellfun(lambdafun,Sigma_x,invB,arrayfun(@(x){x},gamma_old));           
            lambdaComp = sum(cell2mat(lambdatmp));
            LAMBDA     = norm(y - Phi * mu_x, 2)^2/N ...
                         + LAMBDA*(length(mu_x)-lambdaComp)/N;
    end
        

    % ================= Check stopping conditions, eyc. ==============
    if (size(mu_x) == size(mu_old))
        dmu = max(max(abs(mu_old - mu_x)));
        if (dmu < EPSILON), break; end
    end
    if (PRINT) 
        disp([' iters: ',       num2str(count),...
            ' num coeffs: ',    num2str(usedNum), ...
            ' min gamma: ',     num2str(min(gamma)),...
            ' gamma change: ',  num2str(max(abs(gamma - gamma_old))),...
            ' mu change: ',     num2str(dmu)]); 
    end;

end % main while loop

if (count >= MAX_ITERS) && PRINT
    fprintf('Max iterations (%d) reached. BSBL solution may be sub-optimal.\n\n',MAX_ITERS);
end


if isempty(keep_list)
    Result.x            = zeros(M,1);
    Result.gamma_used   = [];
    Result.gamma_est    = zeros(nblks,1);
    Result.B            = B;
    Result.count        = count;
    Result.lambdatrace  = LAMBDA;
else
    %% Expand hyperparameters
    gamma_used          = sort(keep_list);
    gamma_est           = zeros(nblks,1);
    gamma_est(keep_list,1) = gamma;


    %% reconstruct the original signal
    x = zeros(M,1);
    curLoc = 0;
    for i = 1 : usedNum

        curBlkLen   = size(Sigma0{i},1);
        curLoc      = curLoc + 1;
        seg         = curLoc : 1 : curLoc + curBlkLen - 1;

        realLocs    = blkStartLoc0(keep_list(i)) : blkStartLoc0(keep_list(i))+curBlkLen-1;

        x( realLocs ) = mu_x( seg );
        curLoc        = seg(end);
    end

    % restore parameters to original scale
    if (scl < 0.4) || (scl > 1)
        Result.x = x * scl/0.4;
    else
        Result.x = x;
    end
    
    Result.gamma_used   = gamma_used;
    Result.gamma_est    = gamma_est;
    Result.B            = B;
    Result.count        = count;
    Result.lambda       = LAMBDA;
end
