function [J,alpha,beta,T] = dynamicLoreta(UtY,s2,iLV,alpha,beta,L,U)

%[J,varargout] = dynamicLoreta(V,varargin)
%
% Computes the posterior distribution of the parameters J given some data V. 
% The program solves levels of inference: 1) optimization of parameters J, and
% 2) optimization of hyperparameters alpha and beta. See Trujillo-Barreto
% et. al. (2004) for details.
%
% UtY,s2, and iLV are defined as follows: 
%     Y: Nsensors x time points data matrix
%     K: N x P predictor matrix
%     L: sparse P x P square root of the precision matrix 
%     [U,s,V] = svd( K*inv(L) )
%     iLV = inv(L)*V
%     UtY = U'*Y
%     s2  = s.^2
%
% alpha, beta: hyperparameters
% J: estimated parapeters
% 
%                     P(V|J,alpha)*P(J|beta)
% P(J|V,alpha,beta) = ---------------------- 
%                        P(V|alpha,beta)
% 
% Author: Alejandro Ojeda, SCCN/INC/UCSD, Jan-2013
%
% References:
%   Trujillo-Barreto, N., Aubert-Vazquez, E., Valdes-Sosa, P.A., 2004.
%     Bayesian model averaging in EEG/MEG imaging. NeuroImage 21, 1300–1319



if nargin < 3, error('Not enough input arguments.');end

n = size(UtY,1);
p = length(s);
T = [];

% Initialize hyperparameters
if nargin < 5
    nlambda = 100;
    tol = max([n p])*eps(max(s));
    lambda2 = logspace(log10(tol),log10(max(s)),nlambda);
    gcv = zeros(nlambda,1);
    parfor it=1:nlambda
        d = lambda2(it)./(s2+lambda2(it));
        f = diag(d)*UtY;
        gcv(it) = dot(f,f,1)/sum(d)^2;
    end
    loc = getMinima(gcv);
    if isempty(loc), loc = 1;end
    loc = loc(end);
    lambda = lambda2(loc);

    alpha = rand;
    beta = alpha*lambda;
end

% updating parameters
if nargin > 6 && nargout > 3
    T = iLV*spdiags(alpha*s./(alpha*s.^2+beta))*U';
    J = T*Y;
else
    J = iLV*spdiags(alpha*s./(alpha*s.^2+beta))*UtY;
end

% updating hyperparameters
G = L*J;
gamma = p-beta*(alpha*sum(s2)+beta);
alpha = n-gamma;
beta  = gamma/G'*G;
end


%---
function indmin = getMinima(x)
fminor = diff(x)>=0;
fminor = ~fminor(1:end-1, :) & fminor(2:end, :);
fminor = [0; fminor; 0];
indmin = find(fminor);
end