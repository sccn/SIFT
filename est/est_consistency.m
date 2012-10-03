
function PC = est_consistency(data,AR,PE)
%
% Calculate the percent consistency for a fitted VAR[p] model. The
% percent consistency [1] is an index of the ability for a VAR model, fit 
% to data X, to generate data with the same covariance structure as X. 
% If Rr and Rs represent the vectorized corss-correlation matrices of the
% real and simulated data, respectively, then the percent consistency is
% given by:
% 
% PC = (1 - norm(Rs - Rr)/norm(Rr)) * 100
%
% Inputs: 
%
%   data    [nchs x winlen x ntr]
%   AR      MODEL model coefficients as returned by mvar()
%   PE      Noise covariance matrix
%
% Outputs:
%
%   PC      percent consistency [0 100]%
%
% See Also: est_checkMVARConsistency()
%
% References:
% 
% [1] Ding M, Bressler SL, Yang W, Liang H (2000) Short-window spectral 
%   analysis of cortical event-related potentials by adaptive multivariate 
%   autoregressive modeling: data preprocessing, model validation, and 
%   variability assessment. Biol. Cybern. 83:35-45 
% [2] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapters 3,6. 
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% 
% Author: Tim Mullen 2010, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

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


[nchs m] = size(AR);
p = m/nchs;    % model order
if size(PE,2)>nchs
    C = PE(:,fix(nchs*p+1):fix(nchs*(p+1)));
else
    C = PE;
end

nlags = max(20,p);  %p
ndisc = 10^3;   % number of simulated samples to discard
[nchs pnts ntr] = size(data);

% simulate ntr realizations from VAR model
datasim = tvarsim(zeros(1,nchs),AR,C,[pnts*ntr 1],ndisc);
datasim = reshape(datasim',[nchs,pnts,ntr]);
datasim = permute(datasim,[2 1 3]); % [pnts nchs ntr]

% convert to 2-D matrix with p nans between trials
datasim = nanpad(datasim,p);
data    = nanpad(permute(data,[2 1 3]),p);

% pre-normalize data ignoring NaNs and zero out NaNs
datasim = bsxfun(@minus,datasim,nanmean(datasim));
datasim = bsxfun(@rdivide,datasim,nanstd(datasim));
datasim(isnan(datasim)) = 0;

data = bsxfun(@minus,data,nanmean(data));
data = bsxfun(@rdivide,data,nanstd(data));
data(isnan(data)) = 0;

% compute xcorr (taking trial boundaries into account)
Rs = xcorr(datasim,nlags,'coeff');
Rr = xcorr(data,nlags,'coeff');

% % calculate average auto- and cross-correlation coefficients
% Rs=zeros(2*(nlags+1)-1,nchs^2);
% Rr=Rs;
% for tr=1:ntr
%     Rs = Rs+xcorr(datasim(:,:,tr)',nlags,'coeff');
%     Rr = Rr+xcorr(data(:,:,tr)',nlags,'coeff');
% end
% Rs=Rs./ntr;
% Rr=Rr./ntr;

% due to correlation symmetry, we don't need both correlation vectors RiRj and RjRi
nn=nonzeros(tril(reshape(1:nchs*nchs, nchs, nchs))); 
nzlags = [1:nlags nlags+2:2*(nlags+1)-1];
Rs=Rs(nzlags,nn);
Rr=Rr(nzlags,nn);

% Consistency Index
PC = 100 * (1 - (norm(Rs(:)-Rr(:)) / norm(Rr(:))));

