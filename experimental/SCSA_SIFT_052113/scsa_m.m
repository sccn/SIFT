function [H, status] = scsa_m(X, B, x0, lambda, varargin)

[M, N, T] = size(X);
[M, M, K] = size(x0);

S = reshape(B*reshape(X, M, N*T), M, N, T);
Y = vec(reshape(S(:, K+1:end, :), M, [])');

X = [];
for t = 1:K
  X = cat(3, X, reshape(S(:, K+1-t:end-t, :), M, [])'); 
end
X = reshape(permute(X, [1 3 2]), (N-K)*T, M*K);

indsB = vec(repmat((1:K*(M+1):K*M^2), K, 1) + repmat((0:(K-1))', 1, M))';
indsA = setdiff(1:K*M^2, indsB);
[so indsT] = sort([indsA indsB]);

x0 = permute(x0, [3 2 1]);
AA = {@xforth, @xback, (N-K)*T*M, M*M*K};
[H, status] = dalhsgl(vec(x0([indsA indsB])), [], AA, [], Y, lambda, struct('solver', 'cg', 'display', 0, 'blks', [K*ones(M*(M-1), 1); K*M], varargin{:}));

H2 = zeros(M^2*K, 1);
H2(indsA) = vec(H(1:(K*M*(M-1))));
H2(indsB) = vec(H(((K*M*(M-1))+1):end));

H = permute(reshape(full(H2), K, M, M), [3 2 1]);

function xfo = xforth(x)
    xfo = vec(X*reshape(x(indsT), M*K, []));
end
    
function xba = xback(x)
    xba = vec(X'*reshape(x, (N-K)*T, []));
    xba = xba([indsA indsB]);
end

end

