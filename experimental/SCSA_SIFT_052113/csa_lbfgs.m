function [W f S] = csa_lbfgs(X, K, pars)

[M N T] = size(X);
nstart = 1;
x0     = zeros(M*M*(K+1),nstart);
x      = x0;
f      = zeros(1,nstart);

if nargin > 2
  if isfield(pars, 'nstart')
    nstart = pars.nstart;
  end
  if isfield(pars, 'start')
    x0 = reshape(pars.x0, M*M*(K+1), []);
    nstart = size(x0, 2);
  end
end

for istart = 1:nstart
  x0(:, istart) = vec(cat(3, eye(M), zeros(M, M, K)));
  if istart > 1
    x0(:, istart) = x0(:, istart) + 0.1*randn(size(x0(:, istart)));
  end
end

opts = optimset('Display','off', 'GradObj', 'on', 'FinDiffType', 'central', 'DerivativeCheck', 'off');

% tic
for istart = 1:nstart
%   disp(istart)
  [x(:, istart) f(istart)] = fminunc(@(x) csa_fgrad(x), x0(:, istart), opts);
end
% toc

[f in] = min(f);
x = x(:, in);
  
W = reshape(x, M, M, K+1);

if nargout > 2
  S = reshape(W(:, :, 1)*reshape(X, M, []), M, N, T);
end

function [f grad] = csa_fgrad(x)
  
  W = reshape(x, [M M K+1]);

  I = zeros(M, (N-K)*T);
  for k = 0:K
      I = I + W(:, :, k+1)*reshape(X(:, K-k+1:end-k, :), M, []);
  end
  f = (K-N)*log(abs(det(W(:, :, 1)))) - sum(sum(log(sech(I)/pi)));

  if nargout > 1
    grad = zeros(M, M, K+1);
    grad(:, :, 1) = (K-N)*inv(W(:, :, 1))';
    for k = 0:K
        grad(:, :, k+1) = grad(:, :, k+1) + tanh(I)*reshape(X(:, K-k+1:end-k, :), M, [])';
    end
    grad = vec(grad);
  end
end

end