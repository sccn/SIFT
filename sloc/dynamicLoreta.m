function [J,alpha,beta,T] = dynamicLoreta(Ut,Y,s2,iLV,L,options,alpha,beta)

%[J,varargout] = dynamicLoreta(V,varargin)
%
% Computes the posterior distribution of the parameters J given some data V. 
% The program solves levels of inference: 1) optimization of parameters J, and
% 2) optimization of hyperparameters alpha and beta. See Trujillo-Barreto
% et. al. (2004) for details.
%
% Ut,s2, and iLV are defined as follows: 
%     Y: Nsensors x time points data matrix
%     K: N x P predictor matrix
%     L: sparse P x P square root of the precision matrix 
%     [U,s,V] = svd( K*inv(L) )
%     iLV = inv(L)*V
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


if nargin < 5, error('Not enough input arguments.');end
if nargin < 6
    options.maxTol          = 1e-3;
    options.maxIter         = 100;
    options.verbose         = true;
    options.initNoiseFactor = 0.001;
end


s = s2.^(0.5);
n = size(Ut,1);
p = length(s);

% Initialize hyperparameters
if nargin < 7
    UtY     = Ut*Y;
    tol     = max([n p])*eps(max(s));
    nlambda = options.maxIter;
    lambda2 = logspace(log10(tol),log10(max(s)),options.maxIter);
    gcv     = zeros(nlambda,1);
    parfor k=1:nlambda
        d = lambda2(k)./(s2+lambda2(k));
        f = diag(d)*UtY;
        gcv(k) = dot(f,f,1)/sum(d)^2;
    end
    loc = getMinima(gcv);
    if isempty(loc), loc = 1; end
    loc = loc(end);
    lambda2 = lambda2(loc);
     
    alpha = options.initNoiseFactor*(Y(:)'*Y(:))/n;
    beta  = alpha*lambda2;
end
err = inf;

for it=1:options.maxIter
    if err < options.maxTol, break;end
    
    % updating parameters
    % J = iLV*diag(alpha*s./(alpha*s.^2+beta))*UtY;
    T = iLV*diag(alpha*s./(alpha*s2+beta))*Ut;
    J = T*Y;
    
    % updating hyperparameters
    alpha_old = alpha;
    beta_old  = beta;
    gamma     = p-beta*sum(1./(alpha*s2+beta));
    alpha     = n-gamma;
    beta      = gamma/norm(L*J,'fro');
    err       = 0.5*(abs(alpha_old-alpha) + abs(beta_old-beta));
    if options.verbose
        disp([num2str(it) ' => alpha: ' num2str(alpha) '  beta: ' num2str(beta) ' df: ' num2str(gamma) ' error: ' num2str(err)]);
    end
end
if it == options.maxIter, warning('Maximum iterations reached. Failed to converge.');end
end


%---
function indmin = getMinima(x)
fminor = diff(x)>=0;
fminor = ~fminor(1:end-1, :) & fminor(2:end, :);
fminor = [0; fminor; 0];
indmin = find(fminor);
end

% function T = standT(T,gamma,Y,L,iLV,s,Ut)
% H = (L*iLV*diag(s)*Ut)'*T;
% E = Y-H*Y;
% sigma = E'*E/gamma;
% dT = 1./sqrt(dot(T,T,2));
% S = 1./sigma*dT;
% T = bsxfun(@times,S,T);
% end