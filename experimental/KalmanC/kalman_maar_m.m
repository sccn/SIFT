function [x,err,Kalman] = kalman_maar_m(y,p,UC,Kalman,mode)

[tmp,M] = size( y );
L = M*M*p;

if nargin<5
    mode = 0;
end

if nargin<4
    % initialize kalman filter
	Kalman.G	= zeros(L,M);
	Kalman.H	= zeros(M,L);
	Kalman.Kp	= eye( L );
	Kalman.Q1	= eye( L ) * UC;
	Kalman.Q2	= eye(M);
	Kalman.errCov   = eye(M);
	Kalman.x	= zeros(L,1);
	% vector of past observations
	Kalman.yr	= zeros(1,M*p);
    %Kalman.yr = [y(:)' zeros(1,M*(p-1))];
    x = Kalman.x';
    err = zeros(size(y));
    %err = y;
    return
end

% update measurement matrix
%  special case of Kalman.H = kron( eye(M), Kalman.yr );
for m = 1 : M 
    for n = 1 : M*p 
        Kalman.H(m+(n-1)*M+(m-1)*L) = Kalman.yr(n);
    end
end

% calculate prediction error
ye = (Kalman.H*Kalman.x)';
err = y - ye;

if ~any(isnan(err(:)))
	% Adaptive Error covariance
	Kalman.errCov = (1-UC)*Kalman.errCov + UC*(err'*err);

	if( mode == 0 )
		Kalman.Q2 = Kalman.Q2;	% No update
	else
		% update Q2 using the prediction error
		% of the previous step (corresponds to BioSigs's eMode==4)
		Kalman.Q2 = Kalman.errCov;
	end
	
	KpH = Kalman.Kp * Kalman.H';
	HKp = Kalman.H * Kalman.Kp;
	
	% Kalman Gain
	%Kalman.G = KpH * inv( Kalman.H*KpH + Kalman.Q2 );
	Kalman.G = (( Kalman.H*KpH + Kalman.Q2 ) \ KpH')';
	
	% calculation of the a-posteriori state error
	% covariance matrix
	K = Kalman.Kp-Kalman.G*HKp;
	
	% updating Q1
	if( mode == 0 )
		Kalman.Q1 = Kalman.Q1;	% no update
	else
		% corresponds to aMode==17
		K = 0.5 * (K+K');
		Kalman.Q1 = UC * Kalman.Kp;
	end
	
	% a-priori state error covariance matrix
	% for the next time step
	Kalman.Kp=K+Kalman.Q1;
	
	% current estimation of state x
	Kalman.x = Kalman.x + Kalman.G*err';
	
end

% add new observation, drop oldest
Kalman.yr = [y(:)' Kalman.yr(1:M*(p-1))];

x = Kalman.x';

end
