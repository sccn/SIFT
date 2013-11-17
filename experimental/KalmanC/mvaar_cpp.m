function [xout,errCov] = mvaar_cpp(y,p,UC,mode,verb,downsampleFactor,constraints)
% Multivariate (Vector) adaptive AR estimation base on a multidimensional
% Kalman filer algorithm. A standard VAR model (A0=I) is implemented. The
% state vector is defined as X=(A1|A2...|Ap)' and x=vec(X)
%
% [x,e,Kalman,Q2] = mvaar(y,p,UC,mode,Kalman)
%
% The standard MVAR model is defined as:
%
%		y(n)-A1(n)*y(n-1)-...-Ap(n)*y(n-p)=e(n)
%
%	The dimension of y(n) equals M
%
%	Input Parameters:
%
% 		y			Observed data or signal [N x M] where N = epoch length
% 		p			prescribed maximum model order (default 1)
%		UC			update coefficient	(default 0.001)
%		mode	 	update method of the process noise covariance matrix 0...7 ^
%					correspond to S0...S7 (default 0)
%       verb        verbosity
%       downsampleFactor:   Starting from sample t=max(2,k), store only every k Kalman
%                           coefficients (states, etc) where k=downsampleFactor.
%       constraints structure with fields .D and .d containing constraints
%                   of the form Dx = d (see [1] below)
%
%	Output Parameters
%
%		e			prediction error of dimension s
%		x			state matrix of dimension [T x M*M*p]
%                   where T = ceil((N-downsampleFactor+q)/downsampleFactor)
%                   where q = (downsampleFactor>1 ? 1 : 0)
%                   - note that we never store the coefficient matrix for t=1
%                   since this is always zero (we have no sample at t=0 from
%                   which to compute the coefficient matrix for t=1)
%		Q2			measurement noise covariance matrix of dimension M x M
%       Kout        estimated state noise covariance matrix
%       Kalman      Kalman structure (can be used as subsequent startup)
%

%       $Id: mvaar.m 5090 2008-06-05 08:12:04Z schloegl $
%       Copyright (C) 2001-2002 Christian Kasess
%       Copyright (C) 2003, 2008 Alois Schloegl
%
%       Copyright (C) 2010-2011 Tim Mullen
%
%       Modified by Tim Mullen
%       01/23/2011 -- Modified for downsampled storage
%       04/13/2011 -- Optimized performance
%       05/12/2011 -- Added projection onto constraint surface [1]
%       05/20/2011 -- Added additional noise covariance update modes
%
%       [1] Simon D (2010) Kalman Filtering with State Constraints:
%       A Survey of Linear and Nonlinear Algorithms.
%       Control Theory & Applications, IET 4:1303-€“1318
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin<6
    verb = 2;
end
if nargin<4,
    mode=0;
end;
if nargin<3,
    UC=0.001;
end;
if nargin<2,
    p=1;
end
if nargin<1,
    fprintf(2,'No arguments supplied\n');
    return
end;
if nargin<8 || isempty(constraints)
    doConstraints = 0;
else
    doConstraints = 1;
    Constr_D = constraints.D;
    Constr_d = constraints.d;
end

if ~any(mode==(0:8))
    fprintf(2,'Invalid mode (0...8)\n');
    return
end;


[LEN, M, NTR] = size(y);		% signal length, number of channels, number of trials
L = M*M*p;

if LEN<(p+1),
    fprintf(2,'Not enough observed data supplied for given model order\n');
    return
end

% ye = zeros(size(y));	%prediction of y

% size of downsampled storage
dslen = ceil((LEN-downsampleFactor+(downsampleFactor>1))/downsampleFactor);

xout=zeros(L,dslen,NTR);

if nargout>1
    errCov=zeros(M,M,dslen,NTR);
end

if verb==2
    h=waitbar(0,sprintf('fitting VAR[%d] model [mode=%s] ...', ...
        p, 'Kalman'));
end


% initialize Kalman model
[x err Kalman] = kalman_maar(y,p,UC);


for tr=1:NTR
    
    % NOTE: should we re-initialize the state variables here?
    
    curval = 1;
    
    for n = 2:LEN
        
        if verb==2 && ~mod(n,100)
            waitbar(n/LEN,h,...
                {sprintf('Trial (%d/%d)',tr,NTR), ...
                 sprintf('fitting VAR[%d] model [mode=%s] (%d/%d) ...',p,'Kalman',n,LEN)});
        end
        
        if(n<=p)
            Yr=[y(n-1:-1:1,:,tr)' zeros(M,p-n+1)];	%vector of past observations
            Yr=Yr(:)';
        else
            Yr=y(n-1:-1:n-p,:,tr)';				%vector of past observations
            Yr=Yr(:)';
        end
                
        if ~any(isnan(Yr)),
            
            % update model
            [x err Kalman] = kalman_maar(Yr,p,UC,Kalman,mode);
            
            if doConstraints
                KD = K*Constr_D';
                
                % project the solution onto the constraint surface
                x = x - (KD/(Constr_D*KD))*(Constr_D*x - Constr_d);
            end
            
        end; % isnan(err)
        
        if ~mod(n,downsampleFactor)
            % store the current state
            
            if 0; doConstraints;
                
                KD = K*Constr_D';
                
                % project the solution onto the constraint surface
                xout(:,curval,tr) = x - (KD/(Constr_D*KD))*(Constr_D*x - Constr_d);
            else
                
                xout(:,curval,tr) = x;
            end
            
            if nargout>1
                errCov(:,:,curval,tr)=Q2;
            end;
            
            curval = curval + 1;
        end
    end;
    
end

xout = permute(xout,[2 1 3]);

if verb==2, close(h); end

