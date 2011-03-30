
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

% simulate data from model
datasim = zeros(nchs,pnts,ntr);
[datasim]=arsim(zeros(1,nchs),AR,C,pnts*ntr,ndisc);
datasim = reshape(datasim',[nchs pnts ntr]);


% calculate auto- and cross-correlation coefficients
Rs=zeros(2*(nlags+1)-1,nchs^2);
Rr=Rs;
for tr=1:ntr
    Rs = Rs+xcorr(datasim(:,:,tr)',nlags,'coeff');
    Rr = Rr+xcorr(data(:,:,tr)',nlags,'coeff');
end
Rs=Rs./ntr;
Rr=Rr./ntr;

% due to correlation symmetry, we don't need both correlation vectors RiRj and RjRi
nn=nonzeros(tril(reshape(1:nchs*nchs, nchs, nchs))); 
nzlags = [1:nlags nlags+2:2*(nlags+1)-1];
Rs=Rs(nzlags,nn);
Rr=Rr(nzlags,nn);

% Consistency Index
PC = 100 * (1 - (norm(Rs(:)-Rr(:)) / norm(Rr(:))));

