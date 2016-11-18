function [B H f S] = scsa_em(X, K, pars)
%SCSA_EM - sparsely connected sources analysis [1]
%estimates a square demixing B into sources following a sparse MVAR
%model with coefficients H and hyperbolic secant distributed innovations.
%
%The function first learns a non-sparse (unregularized) model via
%maximum-likelihood. This model is then used as a starting point for
%updating B and H in an expectation maximization (EM) fashion, where an
%additional sparsity constraint is imposed on the entries of H.
%
%   inputs : X    - M x N x T data tensor
%            K    - order of source MVAR process
%            pars - optional parameter struct with potential fields
%              .lambda - scalar regularization constant (high values
%                        promote sparsity of the source MVAR model)
%              .emiter - scalar number of EM iterations
%              .nstart - number of initializations
%              .start
%                 .B   - M * M x NSTART tensor of initializations for B
%                 .H   - M x M x K x NSTART tensor of initializations for H
%              .emstop - stopping criterion for E-M (~1e-6 is a good choice).
%                        Algorithm terminates when change in objective func
%                        is less than this
%   outputs: B - M x M demixing matrix
%            H - M x M x K MVAR coefficient tensor
%            f - final value of the SCSA cost function (data likelihood +
%                sparsity penalty)
%            S - ICA source activations
%
%  [1] S. Haufe, R. Tomioka, G. Nolte, K. R. Müller, M. Kawanabe,
%      Modeling sparse connectivity between underlying brain sources
%      for EEG/MEG. IEEE TBME, 57(8):1954 - 1963, 2010
%

[M N T] = size(X);
nstart = 1;
lambda = 0;
emiter = 5;

inddiag = 1:M:M^2;
indoffdiag = setdiff(1:M^2, inddiag);

if nargin > 2
    if isfield(pars, 'runCSA')
        runCSA = pars.runCSA;
    end
    if isfield(pars, 'nstart')
        nstart = pars.nstart;
    end
    if isfield(pars, 'start')
        nstart = size(pars.start.B, 3);
    else
        pars.start.B = repmat(eye(M), [1, 1, nstart]);
        pars.start.H = zeros(M, M, K, nstart);
        pars.start.B(:, :, 2:nstart) = pars.start.B(:, :, 2:nstart) + 0.1*randn(M, M, nstart-1);
        pars.start.H(:, :, :, 2:nstart) = pars.start.H(:, :, :, 2:nstart) + 0.1*randn(M, M, K, nstart-1);
    end
    if isfield(pars, 'lambda')
        lambda = pars.lambda;
    end
    if isfield(pars, 'emiter')
        emiter = pars.emiter;
    end
    if ~isfield(pars,'verb');
        pars.verb = false;
    end
end

% initialize variables
B = zeros([size(pars.start.B), nstart]);
H = zeros([size(pars.start.H), nstart]);
f = zeros(1,nstart);

nops = nstart*emiter;
iop  = 1;
waitbarTitle = 'Running SCSA...';
if pars.verb==1
    disp(waitbarTitle);
elseif pars.verb==2
    multiWaitbar(waitbarTitle, ...
        'Color', [0 0 1], ...
        'CanCancel','on', ...
        'CancelFcn',@(a,b) disp('[Cancel requested. Please wait...]'));
end

qflag = false;
for istart = 1:nstart
    if qflag, break; end
    if pars.verb==1
        fprintf('-------------------------------\n');
        fprintf('ITERATION %d\n',istart);
    end
    if runCSA
        if pars.verb==1
            fprintf('...finding initial solution using CSA...\n'); end
        W = csa_lbfgs(X, K, struct('x0', vec(scsa2csa(pars.start.B(:, :, istart), pars.start.H(:, :, :, istart)))));
        [B(:, :, istart) H(:, :, :, istart)] = csa2scsa(W);
    else
        B(:, :, istart) = pars.start.B(:, :, istart);
        H(:, :, :, istart) = pars.start.H(:, :, :, istart);
    end
    if pars.verb==1
        fprintf('...evaluating cost function...\n'); end
    f(istart) = scsa_eval(B(:, :, istart), H(:, :, :, istart));
    
    if lambda > 0
        for iem = 1:emiter
            if qflag, break; end
            if pars.verb==1
                fprintf('EM Iteration: %0.4f\n',iem/emiter);
            end
            
            f_old = f;
            
            % M-step
            % --------------------------------------------------------------------------------------
            % FIXME: this could be replaced by a call to mvar_dalSCSA
            H(:, :, :, istart) = scsa_m(X, B(:, :, istart), H(:, :, :, istart), lambda);
            
            % E-step
            % --------------------------------------------------------------------------------------
            B(:, :, istart) = scsa_e(X, H(:, :, :, istart), B(:, :, istart));
            f(istart) = scsa_eval(B(:, :, istart), H(:, :, :, istart));
            
            if norm(f_old-f(istart))/norm(f(istart)) < pars.emstop
                if pars.verb
                    fprintf('SCSA early termination: stopping criterion (eps<%0.5g) reached\n',pars.emstop);
                end
                break;
            end
            
            if pars.verb==2
                drawnow;
                % graphical waitbar
                cancel = multiWaitbar(waitbarTitle,iop/nops);
                if cancel
                    if strcmpi('yes',questdlg2( ...
                            'Are you sure you want to cancel?', ...
                            'Model Fitting','Yes','No','No'));
                        multiWaitbar(waitbarTitle,'Close');
                        qflag = true;
                        continue;
                    else
                        multiWaitbar(waitbarTitle,'ResetCancel',true);
                    end
                end
                iop = iop + 1;
            end
        end
    end
    
end


% get parameters that minimize cost function
[f in] = min(f);
B = B(:, :, in);
H = H(:, :, :, in);

if nargout > 3
    S = reshape(B*reshape(X, M, []), M, N, T);
end

if pars.verb==2
    multiWaitbar(waitbarTitle,'Close'); 
end



% Helper Functions
% --------------------------------------------------------------------------------------------------

    function f = scsa_eval(B, H)
        % evaluate the cost function
        S = reshape(B*reshape(X, M, N*T), M, N, T);
        I = reshape(S(:, (K+1):end, :), M, []);
        for k = 1:K
            I = I-H(:, :, k)*reshape(S(:, K-k+1:end-k, :), M, []);
        end
        H = reshape(H, M^2, []);
        
        f = (K-N)*log(abs(det(B))) - sum(sum(log(sech(I)/pi))) ...
            + lambda*(norm(vec(H(indoffdiag, :))) + sum(norms(H(indoffdiag, :), 2, 2)));
    end

end