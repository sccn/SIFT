function Conn = est_mvtransfer(AR,C,freqs,srate,connmethods)
%
% Calculate spectral, coherence, and connectivity measures from a fitted 
% VAR model. See [1] for additional details on VAR model fitting and 
% connectivity measures.  
%
% Input:
%
%   AR:             [M x M*p]   Matrix of autoregressive coefficients for 
%                               M-variate process of order p.
%   C:              [M x M]     Noise covariance matrix
%   freqs:          [vector]    Vector of frequencies of interest (Hz)
%   srate:          [single]    Sampling rate of the data
%   connmethods:    {vector}    Cell vector of strings containing names of
%                               connectivity/spectral estimators. Allowed 
%                               values are:
%           'dDTF':             Direct Directed Transfer Function []
%           'dDTF08':           Modified direct DTF with normalization
%                               across full causal grid
%           'ffDTF':            Full-frequency direct DTF
%               ...
%
% Output:
%
%   Conn:           (struct)    Connectivity structure with following
%                               fields:
%           .<connmethod>       An [M x M x length(freqs) x length(times)]
%                               matrix of connectivity estimates for
%                               one of the estimators named above
%
% Supported measures:
% 
%     'dDTF'
%     'dDTF08'
%     'ffDTF'
%     'nDTF'
%     'GGC'
%     'GGC2'
%     'iCoh'
%     'Coh'
%     'S'
%     'pCoh'
%     'mCoh'
%     'RPDC'
%     'GPDC'
%     'nPDC'
%     'PDCF'
%     'DTF'
%     'Sinv'
%     'PDC'
% 
% See Also: est_mvarConnectivity(), pop_est_mvarConnectivity()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapters 3,6. 
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% 
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu
% This function was inspired by the mvfreqz() function in the TSA toolbox of
% A. Schloegl

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



if ischar(connmethods)
    connmethods = {connmethods};
end

% depedency graph -- each of the measures on the right depend on the
% measure on the left (dependency arrow goes from left to right
% e.g., 'Coh',{'iCoh'} means iCoh depends on Coh (Coh-->iCoh)
% On the RHS, only include measures which are *directly* dependent on 
% the LHS measure. e.g., if A->B->C, then you should have three entries
%  'A', {'B'};
%  'B', {'C'};
%  'C', {}
% If a new measure is added, you MUST update this depedency table
% (add the method descriptor to the LHS, add any children of this method 
% to the RHS and update the RHS of any other entries which this method 
% depends on)
dependencies = {...
                'dDTF',     {};
                'dDTF08',   {};
                'ffDTF',    {'dDTF'};
                'nDTF',     {};
                'GGC',      {};
                'GGC2',     {};
                'iCoh',     {};
                'Coh',      {'iCoh'};
                'S',        {'Coh','GGC','GGC2'};
                'pCoh',     {'dDTF','dDTF08'};
                'mCoh',     {};
                'RPDC',     {};
                'GPDC',     {};
                'nPDC',     {};
                'PDCF',     {};
                'dtf_denom' {'ffDTF','nDTF'};
                'DTF',      {'nDTF','ffDTF','dDTF08','dtf_denom','S'};
                'Sinv',     {'mCoh','pCoh','PDCF'};
                'PDC',      {'DTF','G','nPDC','GPDC','RPDC','PDCF'};
                'Rinv',     {'RPDC'}};
            
% list of all possible intermediate (and final) estimators we might want to
% calculate. 
alldescriptors = dependencies(:,1);

methodsneeded = alldescriptors(find(isneeded(alldescriptors,connmethods,dependencies)));                
univariate_measures = {'mCoh'};  % list of measures that are univariate (nchs x freqs)
singleton_measures = {'DC'};     % list of measures that do not depend on frequency (nchs x 1)

[nchs,ncoeffs] = size(AR);
morder = ncoeffs/nchs;
nfreqs = length(freqs);

z = 2*pi*1i/srate;


I = eye(nchs);
A = [I -AR];
Cinv=inverse(C);

% initialize objects
tmp=setdiff(methodsneeded,univariate_measures);
for i=1:length(tmp)
    % nchs x nchs x freqs (bivariate) measures
    Conn.(tmp{i}) = zeros(nchs,nchs,nfreqs);
end
tmp=intersect(methodsneeded,univariate_measures);
for i=1:length(tmp)
    % nchs x 1 x freqs (univariate) measures
    Conn.(tmp{i}) = zeros(nchs,1,nfreqs);
end
tmp=intersect(methodsneeded,singleton_measures);
for i=1:length(tmp)
    % nchs x 1 measures
    Conn.(tmp{i}) = zeros(nchs,1);
end

% if any(strcmpi('Sinv',methodsneeded))
%     detSinv = zeros(nfreqs,1);
% end
 
ddc2 = diag(diag(C).^(-1/2)); 

for n=1:nfreqs
     
    if any(strcmpi('PDC',methodsneeded))
        % complex non-normalized PDC
        % (Fourier transform of model coefficients)
        for k = 1:morder+1,
            Conn.PDC(:,:,n) = Conn.PDC(:,:,n) + A(:,k*nchs+(1-nchs:0))*exp(-z*(k-1)*freqs(n));
        end
    end
    
%     imagesc(abs(Conn.PDC(:,:,n).*exp(1i*z)));
    
    if any(strcmpi('DTF',methodsneeded))
        % complex non-normalized DTF
        Conn.DTF(:,:,n)  = inverse(Conn.PDC(:,:,n));
    end
    
    % --- SPECTRAL MEASURES ---
    
    if any(strcmpi('S',methodsneeded))
        % complex spectral matrix
        Conn.S(:,:,n)  = Conn.DTF(:,:,n)*C*Conn.DTF(:,:,n)';  %/srate
        %         S1(:,:,n) = h(:,:,n)*h(:,:,n)';
    end
    
    if any(strcmpi('Sinv',methodsneeded))
        % inverse spectral matrix:  
        % inv(S) = inv(DTF*C*DTF') = inv(DTF)'inv(C)inv(DTF) = PDC'inv(C)PDC
    	Conn.Sinv(:,:,n) = Conn.PDC(:,:,n)'*Cinv*Conn.PDC(:,:,n);
%         detSinv(n) = det(Conn.Sinv(:,:,n));    
    end
    
    if any(strcmpi('Coh',methodsneeded))
        % complex coherency Cxy = Sxy/sqrt(abs(Sxx*Syy))
        autospect = diag(Conn.S(:,:,n));
        Conn.Coh(:,:,n) = Conn.S(:,:,n)./sqrt(repmat(autospect,[1 nchs]).*repmat(autospect',[nchs 1]));
    end
    
    if any(strcmpi('pCoh',methodsneeded))
        % complex partial coherency
        autospect = diag(Conn.Sinv(:,:,n));
        Conn.pCoh(:,:,n) = Conn.Sinv(:,:,n)./sqrt(repmat(autospect,[1 nchs]).*repmat(autospect',[nchs 1]));
    end
    
    if any(strcmpi('iCoh',methodsneeded))
        Conn.iCoh(:,:,n) = imag(Conn.Coh(:,:,n));
    end
    % --- PARTIAL DIRECTED COHERENCE MEASURES ---
    
    if any(strcmpi('GPDC',methodsneeded))
        % generalized PDC
        gtmp = abs(ddc2*Conn.PDC(:,:,n));
        gtmp_denom =diag(gtmp'*gtmp)';
        Conn.GPDC(:,:,n) = gtmp./sqrt(gtmp_denom(ones(1,nchs),:));
    end
    
    if any(strcmpi('nPDC',methodsneeded))
        % complex normalized PDC (measures interactions with respect to
        % given signal source -- worse for identifying sources)
        pdc_denom = diag(abs(Conn.PDC(:,:,n))'*abs(Conn.PDC(:,:,n)))';
        Conn.nPDC(:,:,n) = Conn.PDC(:,:,n)./sqrt(pdc_denom(ones(1,nchs),:));
    end
        
    if any(strcmpi('PDCF',methodsneeded))
        % partial directed coherence factor
        pdcf_denom = diag(Conn.Sinv(:,:,n)).';
        Conn.PDCF(:,:,n) = abs(Conn.PDC(:,:,n))./sqrt(pdcf_denom(ones(1,nchs),:));
    end
        
    % --- DIRECTED TRANSFER FUNCTION MEASURES ---
    
    if any(strcmpi('dtf_denom',methodsneeded))
        tmp = diag(abs(Conn.DTF(:,:,n))*abs(Conn.DTF(:,:,n))');
        Conn.dtf_denom(:,:,n) = tmp(:,ones(1,nchs));  % equivalent to repmat
    end
    
    if any(strcmpi('nDTF',methodsneeded))
        % complex normalized DTF (measures interactions with respect to
        % given signal sink -- worse for identifying sinks)
        Conn.nDTF(:,:,n) = Conn.DTF(:,:,n)./sqrt(Conn.dtf_denom(:,:,n));
    end
    
%     % store the estimates for this frequency
%     for m=methodsneeded
%         eval(sprintf('Conn.(%s)(:,:,n) = %s;',m{1},m{1}));
%     end
    
 end  % for each frequency
 
 
 % --- MORE DTF MEASURES ---
if any(strcmpi('ffDTF',methodsneeded))
    % (real) full-frequency DTF
    Conn.ffDTF = abs(Conn.DTF)./abs(repmat(sqrt(sum(Conn.dtf_denom,3)),[1 1 nfreqs]));
end
 
if any(strcmpi('dDTF',methodsneeded))
    % (real) direct-DTF
    Conn.dDTF = abs(Conn.pCoh).*Conn.ffDTF;
end

if any(strcmpi('dDTF08',methodsneeded))
    % complex direct-DTF with new normalization
    % normalizing by entire causal map
    % (Korziniewska, ..., Crone, et al 2008)
    Conn.dDTF08 = abs(Conn.DTF).*abs(Conn.pCoh);
    Conn.dDTF08 = Conn.dDTF08./sqrt(sum(Conn.dDTF08(:).^2));
end

% --- MORE PDC MEASURES ---
if any(strcmpi('Rinv',methodsneeded))
    % inverse covariance matrix
    Conn.Rinv = est_calcInvCovMat(AR,C);
end

if any(strcmpi('RPDC',methodsneeded))
    V = est_calcInvCovMatFourier(Conn.Rinv,C,freqs,srate,nchs,morder,0);
    for i=1:nchs
        for j=i:nchs
            aij = [real(squeeze(Conn.PDC(i,j,:))).'; imag(squeeze(Conn.PDC(i,j,:))).'];
            aji = [real(squeeze(Conn.PDC(j,i,:))).'; imag(squeeze(Conn.PDC(j,i,:))).'];
            for k=1:nfreqs
                Conn.RPDC(i,j,k) = ((aij(:,k)'/squeeze(V(i,j,k,:,:))))*aij(:,k);
                Conn.RPDC(j,i,k) = ((aji(:,k)'/squeeze(V(j,i,k,:,:))))*aji(:,k);
%                 %% DEBUG!
%                  try chol(squeeze(V(j,i,k,:,:))); % check for pos-definiteness
%                  catch, fprintf('V<0 i=%d j=%d f=%d || ',i,j,k); end
%                  if i~=j && RPDC(i,j,k) < 0 || RPDC(j,i,k) < 0, keyboard; end
%                 %% DEBUG!
            end
        end
    end
end

if any(ismember({'GGC','GGC2'},methodsneeded))
    absS = abs(Conn.S);
    absHsq = abs(Conn.DTF).^2;
end

if any(strcmpi('GGC',methodsneeded))
    % Bivariate Granger-Geweke Causality (similar to Kaminski et al. 2001. )
    for i=1:nchs
        for j=1:nchs
            Conn.GGC(i,j,:) = ((C(i,i)*C(j,j)-C(i,j)^2)/C(j,j))*absHsq(i,j,:)./absS(j,j,:);
        end;
    end;
end

if any(strcmpi('GGC2',methodsneeded))
    % Bivariate Granger-Geweke Causality (similar to Bessler et al. 2007)
    for i=1:nchs
        for j=1:nchs
            Conn.GGC2(i,j,:) = log( absS(i,i,:)./ ...
                                 (absS(i,i,:) - (C(j,j)-(C(i,j)^2/C(i,i)))*absHsq(i,j,:)));
        end;
    end;
end

if any(strcmpi('mCoh',methodsneeded))
    % complex multiple coherence
    detS = zeros(nfreqs,1);
    for f=1:nfreqs, detS(f) = det(Conn.S(:,:,f)); end
    for i=1:nchs
        Conn.mCoh(i,:) = sqrt(ones(nfreqs,1)-(detS./(squeeze(Conn.S(i,i,:)).*hlp_minor(Conn.S,i,i))));
    end
end


% return only the measures requested
Conn = rmfield(Conn,setdiff(fieldnames(Conn),connmethods));

% use single-precision to save space (the double() is to force inverse
% objects to be converted to numerics)
Conn = structfun(@(x) (single(double(x))),Conn,'uniformoutput',false);


%% HELPER FUNCTIONS 

 function needed = isneeded(curmethod,allmethods,dependencies)
% check whether a particular measure will need to be calculated
% Inputs:
%   - curmethod:  a string or cell array of strings containing method
%                 descriptors to check whether needed
%   - allmethods: a cell array of strings containing method descriptors we
%                 want to calculate
% Output:
%   - needed:     array of same size as curmethod containing 1 if method is 
%                 needed, 0 otherwise

nmethods = length(curmethod);

if iscell(curmethod) && nmethods >1
    needed = zeros(1,nmethods);
    for i=1:nmethods
        needed(i) = isneeded(curmethod{i},allmethods,dependencies);
    end
    return;
end
             
if isempty(curmethod)
    needed =0;
    return;
elseif ismember(curmethod,allmethods)
    needed = 1;
    return;
else
    % recursively determine whether any of the other methods the user 
    % wants depend on curmethod
    idx=strmatch(curmethod,dependencies(:,1),'exact');
    if isempty(idx) || isempty(dependencies(idx,2)), needed=0; return; end
    children=dependencies(idx,2);
    if any(isneeded(children{1},allmethods,dependencies))
        needed = 1;
        return;
    end
end

% if we get here, this particular curmethod is not needed
needed = 0;
    
    