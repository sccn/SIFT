function [J,lambdaOpt,Yhat] = ridgeGCV(Y,K,L,lambdaOpt,nlambda,plotGCV,verb,lambdaLearnRule)
%[J,lambdaOpt] = ridgeGCV(Y,K,L,nlambda,plotGCV)
%
% Estimates a ridge regression model, also know as Tikhonov regularization, 
% or minimum norm with L2 prior (or Loreta in the EEG inverse solution literature). 
% For an implementation of sLORETA model see the function inverseSolutionLoreta.
%
% Y: measurements (Nsensors X 1)
% K: N X P predictor matrix
% L: P X P prior covariance matrix (sparse matrix is recomended)
% Jest: estimated parapeters
% nlambda: maximum size of the grid for the hyperparameter lambda, default: 100
% plotGCV: plot the GCV curve (true/false), default: false
% 
% Jest = argmin(J) ||Y-K*J||^2 + lambda*||L*J||^2
% with lambda > 0
%
% This code is based on a previous implementation used in Valdes-Hernandez 
% et al. (2009), written by Alejandro Ojeda and Pedro Valdez-Hernandez at 
% the Cuban Neuroscience Center in 2009.
% 
% Author: Alejandro Ojeda, SCCN/INC/UCSD, Jul-2012
%
% References:
%   Pedro A. Valdés-Hernández, Alejandro Ojeda, Eduardo Martínez-Montes, Agustín
%       Lage-Castellanos, Trinidad Virués-Alba, Lourdes Valdés-Urrutia, Pedro A.
%       Valdes-Sosa, 2009. White matter architecture rather than 
%       cortical surface area correlates with the EEG alpha rhythm. NeuroImage 49
%       (2010) 2328–2339

if nargin < 2, error('Not enough input arguments.');end
[n,p] = size(K);
if nargin < 3 || isempty(L), L = speye(p);end
if nargin < 4 || isempty(lambdaOpt), lambdaOpt = -1; end
if nargin < 5, nlambda = 100;end
if nargin < 6, plotGCV = false;end
if nargin < 7, verb = 0; end
if nargin < 8, lambdaLearnRule = 'grid_gcv'; end
    
[U,S,V] = svd(K/L,'econ');
V = L\V;
s = diag(S);
s2 = s.^2;
UtY = U'*Y;

if lambdaOpt<0 && strcmpi(lambdaLearnRule,'grid_gcv')
    % automatically find optimal lambda via grid search
    % for minimum of GCV
    tol = max([n p])*eps(max(s));
    lambda2 = logspace(log10(tol),log10(max(s)),nlambda);
    gcv = zeros(nlambda,1);
    parfor it=1:nlambda
        d = lambda2(it)./(s2+lambda2(it));
        f = diag(d)*UtY;
        gcv(it) = dot(f,f,1)/sum(d)^2;
    end
    loc = getMinima(gcv);
    if isempty(loc), 
        % no minimum found, search for elbow instead
        if verb
            fprintf('no GCV minimum, finding elbow...\n');
        end
        [dummy loc] = min(gcv); %hlp_findElbow(gcv);
    end
    loc = loc(end);
    lambdaOpt = lambda2(loc);
    if verb
        fprintf('lambdaOpt: %0.5g, GCV: %0.5g\n',lambdaOpt,gcv(loc));
    end
end

T = V*diag(s./(s2+lambdaOpt))*U';
J = T*Y;                            % J = (K'*K+lambda*L'*L)\K'*Y
if nargout > 2, Yhat = K*J;end
if plotGCV
    figure;
    semilogx(lambda2,gcv)
    %plot(lambda2,gcv)
    xlabel('log-lambda');
    ylabel('GCV');
    hold on;
    plot(lambdaOpt,gcv(loc),'rx','linewidth',2)
    grid on;
end

%---
function indmin = getMinima(x)
fminor = diff(x)>=0;
fminor = ~fminor(1:end-1, :) & fminor(2:end, :);
fminor = [0; fminor; 0];
indmin = find(fminor);