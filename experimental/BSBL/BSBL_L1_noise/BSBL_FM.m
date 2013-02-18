function Result = BSBL_FM(PHI,y,blkStartLoc,LearnLambda,varargin)
%----------------------------------------------------------------
% Input for BSBL-FM:
%   PHI: projection matrix
%   y:   CS measurements
%   blkStartLoc : Start location of each block
%   LearnLambda : (1) If LearnLambda = 1, 
%                     use the lambda learning rule for MEDIUM SNR cases (SNR<=30dB)
%                     (using lambda=std(y)*1e-2 or user-input value as initial value)
%                 (2) If LearnLambda = 2, 
%                     use the lambda learning rule for HIGH SNR cases (SNR>30dB) 
%                     (using lambda=std(y)*1e-3 or user-input value as initial value)
%                 (3) If LearnLambda = 0, do not use the lambda learning rule 
%                     ((using lambda=std(y)*1e-5 or user-input value as initial value)
%
% [varargin values -- in most cases you can use the default values]
%   'LEARNTYPE'  : LEARNTYPE = 0: Ignore intra-block correlation
%                  LEARNTYPE = 1: Exploit intra-block correlation 
%                 [ Default: LEARNTYPE = 1 ]
%   'VERBOSE'    : debuging information.
%   'r'          : user defined Toeplitz Coeff
%   'EPSILON'    : convergence criterion
%
% ==============================  OUTPUTS ============================== 
%   Result : 
%      Result.x          : the estimated block sparse signal
%      Result.gamma_used : indexes of nonzero groups in the sparse signal
%      Result.gamma_est  : the gamma values of all the groups of the signal
%      Result.B          : the final mean value of each correlation block
%      Result.count      : iteration times
%      Result.lambda     : the final value of lambda
%
% Coded by Benyuan Liu, 2012
% 
% Reference:
%   [1] Benyuan Liu, Zhilin Zhang, Hongqi Fan, Zaiqi Lu, Qiang Fu, The Fast
%   Marginalized Block SBL Algorithm, submitted to IEEE Signal Processing
%   Letters, 2012
%   [2] Zhilin Zhang, Bhaskar D. Rao, Extension of SBL Algorithms for the 
%   Recovery of Block Sparse Signals with Intra-Block Correlation, to
%   appear in IEEE Transaction on Signal Processing
%
%


% default values for BSBL-FM
eta = 1e-4;      % default convergence test
verbose = 0;     % print some debug information
learnType = 1;   % default to exploit intra block correlation
r = 0;           % default Toeplitz parameter for B and U
if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'learntype'
				learnType = varargin{i+1};
            case 'verbose'    
                verbose = varargin{i+1}; 
            case 'rthd'    
                rthd = varargin{i+1}; 
			case 'epsilon'
				eta = varargin{i+1};
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. intialize, scale
scl = std(y);
if (scl < 0.4) || (scl > 1)
    y = y/scl*0.4;
end
[N,M] = size(PHI);

% select sigma2
sigma2 = 0.01*scl;
if LearnLambda == 0 
	sigma2 = 1e-5*scl;   % noiseless
elseif LearnLambda == 2
	sigma2 = 1e-3*scl;   % high SNR
elseif LearnLambda == 1
	sigma2 = 1e-2*scl;   % medium SNR
end
beta = 1/sigma2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. formalize the blocks and quantities used in the code
%    p           : the number of blocks
%    blkStartLoc : the start index of blk
%    blkLenList  : the length of each block
p = length(blkStartLoc);
blkLenList = ones(p,1);
for k = 1 : p-1
    blkLenList(k) = blkStartLoc(k+1)-blkStartLoc(k);
end
blkLenList(p) = M - blkStartLoc(end)+1;
% 2. prepare the quantities used in the code.
for k = 1 : p
	currentLoc = blkStartLoc(k);
	currentLen = blkLenList(k);
	currentSeg = currentLoc:1:currentLoc + currentLen - 1;

	Phi_k = PHI(:,currentSeg);
	S{k} = beta.*Phi_k'*Phi_k;
	Q{k} = beta.*Phi_k'*y;
end
% 3. start from *NULL*, decide which one to add ->
theta = zeros(p,1);
for k = 1 : p
	A{k} = inv(S{k})*(Q{k}*Q{k}' - S{k})*inv(S{k});
	theta(k) = 1/blkLenList(k) * trace(A{k});
	A{k} = eye(blkLenList(k)).*theta(k);
end
% select the basis that minimize the change of *likelihood* ->
ml = inf*ones(1,p);
ig0 = find(theta>0);
len = length(ig0);
for kk = 1:len
	k = ig0(kk);
	ml(k) = log(abs(det(eye(blkLenList(k)) + A{k}*S{k}))) ...
	      - Q{k}'/(eye(blkLenList(k)) + A{k}*S{k})*A{k}*Q{k};
end
[foo,index] = min(ml);
gamma = theta(index);
Am{index} = A{index}; % Am -> record the past value of A
if verbose, fprintf(1,'ADD,\t idx=%3d, GAMMA_OP=%f\n',index,gamma); end
% 3. update quantities (Sig,Mu,S,Q,Phiu) 
Sigma_ii = (eye(blkLenList(index))/Am{index} + S{index})\eye(blkLenList(index));
Sig = Sigma_ii;
Mu = Sigma_ii*Q{index};
% The relevent block basis
currentLoc = blkStartLoc(index);
currentLen = blkLenList(index);
currentSeg = currentLoc:1:currentLoc + currentLen - 1;
Phi_i = PHI(:,currentSeg); 
Phiu = Phi_i;
for k = 1 : p
	currentLoc = blkStartLoc(k);
	currentLen = blkLenList(k);
	currentSeg = currentLoc:1:currentLoc + currentLen - 1;

	Phi_k = PHI(:,currentSeg);
	S{k} = S{k} - beta^2.*Phi_k'*Phi_i*Sigma_ii*Phi_i'*Phi_k;
	Q{k} = Q{k} - beta.*Phi_k'*Phi_i*Mu;
end
%
max_it = 1000; 
ML=zeros(max_it,1);
r = zeros(p,1);      % record the value of correlations
for count = 1:max_it
	% if we want to learn the intra-block-correlation
	if learnType~=0
		len = length(index);
		for i = 1 : len
			k = index(i);  % k is GLOBAL, i is LOCAL.
			localSeg = extractLocal(index,index(i),blkLenList);
			Sigma_ii = Sig(localSeg,localSeg);
			Mu_i = Mu(localSeg);
			[foo,bar,r(k)] = learnB(Sigma_ii,Mu_i,gamma(i),1);
		end
		r_hat = mean(r(index));
	end
	% calculate s,q
	for k = 1 : p
		which = find(index==k,1);
        if isempty(which)
			s{k} = S{k};
			q{k} = Q{k};
        else
			invDenom = (eye(blkLenList(k)) - S{k}*Am{k})\eye(blkLenList(k));
			s{k} = invDenom*S{k};
			q{k} = invDenom*Q{k};
        end
		% correlation(1) or not(0)
		if learnType == 0
			A{k} = inv(s{k})*(q{k}*q{k}' - s{k})*inv(s{k});
 			theta(k) = 1/blkLenList(k) * trace(A{k});
            % r(k) = 1.0;
			A{k} = eye(blkLenList(k))*theta(k);
        else
			A{k} = inv(s{k})*(q{k}*q{k}' - s{k})*inv(s{k});
			theta(k) = 1/blkLenList(k) * trace(A{k});
			% rr = mean(diag(A{k},1))/mean(diag(A{k}));
			% if abs(rr)>0.99, rr = 0.99*sign(rr); end
			% r(k) = rr;
			Bc = genB(r_hat,blkLenList(k));
			A{k} = Bc.*theta(k);
		end
	end

    % choice the next basis that [minimizes] the cost function
    ml =  inf*ones(1,p);
    ig0 = find(theta>0);
    % index for re-estimate
    [ire,foo,which] = intersect(ig0,index);
    if ~isempty(ire)
		len = length(which);
		for kk = 1:len
			k = ire(kk); m = which(kk);
			ml(k) = log(abs(det(eye(blkLenList(k)) + A{k}*s{k}))) ...
			        -q{k}'/(eye(blkLenList(k)) + A{k}*s{k})*A{k}*q{k} ...
				   -(log(abs(det(eye(blkLenList(k))+ Am{k}*s{k}))) ...
				    -q{k}'/(eye(blkLenList(k)) + Am{k}*s{k})*Am{k}*q{k});
		end
    end
    % index for adding
    iad = setdiff(ig0,ire);
    if ~isempty(iad)
		len = length(iad);
		for kk = 1:len
			k = iad(kk);
			ml(k) = log(abs(det(eye(blkLenList(k)) + A{k}*s{k}))) ...
				   -q{k}'/(eye(blkLenList(k)) + A{k}*s{k})*A{k}*q{k};
		end
    end
    % index for deleting
    is0 = setdiff([1:p],ig0);
    [ide,foo,which] = intersect(is0,index);
    if ~isempty(ide)
		len = length(which);
		for kk = 1:len
			k = ide(kk); m = which(kk);
			ml(k) = -(log(abs(det(eye(blkLenList(k)) + Am{k}*s{k}))) ...
			         -q{k}'/(eye(blkLenList(k)) + Am{k}*s{k})*Am{k}*q{k});
		end
    end

	% as we are minimizing the cost function :
    [ML(count),idx] = min(ml);
    
    % check if terminates?
	if ML(count)>=0, break; end
    if count > 2 && abs(ML(count)-ML(count-1)) < abs(ML(count)-ML(1))*eta, break; end

    % update block gammas
    which = find(index==idx);
	% processing the quantities update
	if ~isempty(which)  % the select basis is already in the *LIST*
		localSeg = extractLocal(index,idx,blkLenList);
		Sig_j = Sig(:,localSeg);
		Sig_jj = Sig(localSeg,localSeg);
		if theta(idx)>0
			%%%% re-estimate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if verbose,fprintf(1,'REE,\t idx=%3d, GAMMA_OP=%f\n',idx,theta(idx));end
			gamma_new = theta(idx);
			ki = Sig_j/(Sig_jj + Am{idx}/(Am{idx} - A{idx})*A{idx})*Sig_j';
			Sig = Sig - ki;
			Mu = Mu - beta.*ki*Phiu'*y;
            for k = 1:p
				currentLoc = blkStartLoc(k);
				currentLen = blkLenList(k);
				currentSeg = currentLoc:1:currentLoc + currentLen - 1;

				Phi_m = PHI(:,currentSeg);
				S{k} = S{k} + beta^2.*Phi_m'*Phiu*ki*Phiu'*Phi_m;
				Q{k} = Q{k} + beta^2.*Phi_m'*Phiu*ki*Phiu'*y;
            end
			%
			gamma(which) = gamma_new; % 1
			Am{idx} = A{idx};         % 2
		else 
			%%%% delete %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if verbose,fprintf(1,'DEL,\t idx=%3d, GAMMA_OP=%f\n',idx,gamma(which));end
			ki = Sig_j/Sig_jj*Sig_j';
			Sig = Sig - ki;
			Mu = Mu - beta.*ki*Phiu'*y;
			for k = 1:p
				currentLoc = blkStartLoc(k);
				currentLen = blkLenList(k);
				currentSeg = currentLoc:1:currentLoc + currentLen - 1;

				Phi_m = PHI(:,currentSeg);
				S{k} = S{k} + beta^2.*Phi_m'*Phiu*ki*Phiu'*Phi_m;
				Q{k} = Q{k} + beta^2.*Phi_m'*Phiu*ki*Phiu'*y;
			end
			% delete relevant basis and block
            index(which) = [];
			Mu(localSeg) = [];
			Sig(:,localSeg) = [];
			Sig(localSeg,:) = [];
			Phiu(:,localSeg) = [];
			%
			gamma(which) = [];       % 1
			Am{idx} = [];            % 2
			invAm{idx} = [];         % 3
		end
	else
		if theta(idx)>0
			%%%% add %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if verbose,fprintf(1,'ADD,\t idx=%3d, GAMMA_OP=%f\n',idx,theta(idx));end
			gamma_new = theta(idx);
			currentLoc = blkStartLoc(idx);
			currentLen = blkLenList(idx);
			currentSeg = currentLoc:1:currentLoc + currentLen - 1;
			Phi_j = PHI(:,currentSeg);
			%
			Sigma_ii = (eye(blkLenList(idx))+A{idx}*S{idx})\A{idx};
			mu_i = Sigma_ii*Q{idx};
			Sigma_11 = Sig + beta^2.*Sig*Phiu'*Phi_j*Sigma_ii*Phi_j'*Phiu*Sig;
			Sigma_12 = -beta.*Sig*Phiu'*Phi_j*Sigma_ii;
			Sigma_21 = Sigma_12';
			mu_1 = Mu - beta.*Sig*Phiu'*Phi_j*mu_i;
			e_i = Phi_j - beta.*Phiu*Sig*Phiu'*Phi_j;
			for k = 1:p
				currentLoc = blkStartLoc(k);
				currentLen = blkLenList(k);
				currentSeg = currentLoc:1:currentLoc + currentLen - 1;

				Phi_m = PHI(:,currentSeg);
				S{k} = S{k} - beta^2.*Phi_m'*e_i*Sigma_ii*e_i'*Phi_m;
				Q{k} = Q{k} - beta.*Phi_m'*e_i*mu_i;
			end
			% adding relevant basis
			Sig = [Sigma_11 Sigma_12; ...
			       Sigma_21 Sigma_ii];
			Mu  = [mu_1; ...
			       mu_i];
			Phiu = [Phiu Phi_j];
			index = [index;idx];
			%
			gamma = [gamma;gamma_new];  % 1
			Am{idx} = A{idx};           % 2
		end
	end

end
if verbose, fprintf(1,'\n'); end;
% format the output ===> X the signal
weights = zeros(M,1);
len = length(index);
formatSeg = [];
for i = 1:len
	k = index(i);
	currentLoc = blkStartLoc(k);
	currentLen = blkLenList(k);
	currentSeg = currentLoc:1:currentLoc + currentLen - 1;

	formatSeg = [formatSeg currentSeg];
end
weights(formatSeg) = Mu;
if (scl < 0.4) || (scl > 1)
    Result.x = weights * scl/0.4;
else
    Result.x = weights;
end
if learnType == 0, Result.r = 0; else, Result.r = r_hat; end
Result.gamma_used = index;
Result.gamma_est = gamma;
Result.count = count;
Result.lambda = sigma2;
% END %

%% sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions of estimating the AR(1) coefficient r and 
% reconstruction the covariance matrix with B^{-1} valid
function [B,U,r] = learnB(Sig,Mu,gamma,type)
len = length(Mu);
if type == 0
	B = eye(len);
	U = eye(len);
	r = 0;
else
	B = (Sig + Mu*Mu')./gamma;
	r = (mean(diag(B,1))/mean(diag(B)));
	if abs(r) >= 0.99, r = 0.99*sign(r); end;
	[B,U] = genB(r,len);
end
% generate B,U according to r,len 
% NOTE: abs(r) should be less than 0.99
function [B,U] = genB(r,len)
bs = [];
for j = 1:len
	bs = [bs r^(j-1)];
end
B = toeplitz(bs);
U = chol(B);

%% sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract local segment within currently included basis
% input:
%      index : current include basis
%      idx   : which basis to be extract. (must be within index)
%      blkLenList.
% output:
%      localSeg : range of basis correspond to idx block in index
function localSeg = extractLocal(index,idx,blkLenList)
len = length(index);
incrementalLoc = 1;
for i = 1:len
	if index(i)==idx
		localLoc = incrementalLoc;
		localLen = blkLenList(idx);
		localSeg = localLoc:1:localLoc + localLen -1;
		return;
	end
	incrementalLoc = incrementalLoc + blkLenList(index(i));
end

