function [B, f] = scsa_e(X, H, x0)

[M, N, T] = size(X);
K = size(H, 3);

if nargin < 3
  x0 = eye(M);
end

opts = optimset('Display','off', 'GradObj', 'on', 'FinDiffType', 'central', 'DerivativeCheck', 'off');

[x f] = fminunc(@(x) scsa_e_fgrad(x), vec(x0), opts);

B = reshape(x, M, M);

function [f grad] = scsa_e_fgrad(x)
  
  B = reshape(x, [M, M]);
  S = reshape(B*reshape(X, M, N*T), M, N, T);
  I = reshape(S(:, (K+1):end, :), M, []);
  for k = 1:K
    I = I-H(:, :, k)*reshape(S(:, K-k+1:end-k, :), M, []);
  end

  f = (K-N)*log(abs(det(B))) - sum(sum(log(sech(I)/pi)));
  
  if nargout > 1
    grad = (K-N)*inv(B)' + tanh(I)*reshape(X(:, (K+1):end, :), M, [])';
    %     
    %     for t = 1:N-K
    %        for m = 1:M
    %           for k = 1:K
    %             grad = grad - (tanh(I(m, t)) * (H(m, :, k)'*X(:, t-k+K)'));
    %           end
    %        end
    %     end
    for m = 1:M
      for k = 1:K
        grad = grad - tprod(tanh(I(m, :))', [-1], tprod(H(m, :, k)', [1], reshape(X(:, (1:N-K)-k+K, :), M, []), [2 3]), [1 2 -1]);
      end
    end
    grad = vec(grad);
  end
end

end

