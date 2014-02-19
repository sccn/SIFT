function [data, Arsig, x, lambdamax]=gen_ar_sech(M, N, P, K, perc)

% P=10; %order of AR-model
% N=1000; %number of data-points
% M=10; %number of channels;
sigma=.2; %scale of random AR-parameters
% K = ceil(M^2/10);

ndisc=1000; %length of ignored start 

if nargin < 5
  perc = 1;
end

% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

if length(K) == 1
    inddiag = linspace(1, M^2, M);
    indndiag = setdiff(1:M^2, inddiag);
    per = randperm(length(indndiag));
    indndiag = indndiag(per);
    indndiag = indndiag(1:K);
    ind = [inddiag indndiag];
else
    inddiag = linspace(1, M^2, M);
    ind = unique([inddiag K]);
end

lambdamax=10;
while lambdamax> 1
    Arsig=[];
    for k=1:P
      aloc = zeros(M);
      aloc(ind) = double(rand(length(ind), 1) < perc).*randn(length(ind), 1)*sigma;
      Arsig=[Arsig,aloc];
    end
    E=eye(M*P);AA=[Arsig;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
%     disp(lambdamax)
    seed = round(sum(100*clock));
    for ij=1:M;
        for jj = 1:(N+ndisc)
            [x(ij, jj), seed] = sech_sample(0, 1, seed); 
        end
    end
    y=x;
    for i=P+1:N+ndisc;
        yloc=reshape(fliplr(y(:,i-P:i-1)),[],1);
        y(:,i)=Arsig*yloc+x(:,i);
    end
    data=y(:,ndisc+1:end);
end

x = x(:, ndisc+1:end);

Arsig = reshape(Arsig, M, M, P);




