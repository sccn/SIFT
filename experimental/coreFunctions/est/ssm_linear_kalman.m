function [VAR,PE,residuals,Kalman,Kout] = ssm_linear_kalman(varargin)
% low-level function for VAR estimation using Kalman filtering
%
%
% Based on the function mvaar.m by Christian Kasess and Alois Schloegl
% (released under GPL in the Biosig toolbox).
%
% Revision history:
%
%       Initial function: $Id: mvaar.m 5090 2008-06-05 08:12:04Z schloegl $
%       Copyright (C) 2001-2002 Christian Kasess
%       Copyright (C) 2003, 2008 Alois Schloegl
%       Copyright (C) 2010-2011 Tim Mullen
%           01/23/2011 -- Modified for downsampled storage
%           04/13/2011 -- Optimized performance
%           05/12/2011 -- Added projection onto constraint surface [1]
%           05/20/2011 -- Added additional noise covariance update mode
%           27/06/2012 -- New input structure and function name
% References:
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

g = arg_define([0 Inf],varargin, ...
        arg_norep({'y','Data','data'},mandatory,[],'[time x chans x trials] data matrix'), ...
        arg({'p','morder','ModelOrder','p'},10,[],'VAR model order'), ...
        arg({'UC','UpdateCoefficient','updatecoeff'},0.001,[0 1],...
            ['Update coefficient.' 'This represents how strongly we weight more recent samples versus older samples when updating the VAR model coefficients. This parameter is inversely proportional to the memory length via the equation:' ... 
            sprintf('\n') ...
            'MemoryLengthSamples = -1/log(1-UpdateCoefficient)' ...
            sprintf('\n') ...
            'Increasing UpdateCoefficient decreases the "memory" of the filter allowing the filter to more rapidly adapt to non-stationary fluctuations in the process mean. This generally comes at a cost of increased variability and decreased stability.']), ...
        arg({'mode','UpdateMode','updatemode'},1,{1 2 3 4 5},'Noise covariance matrix update mode'), ...
        arg_nogui({'Kalman','KalmanObject'},struct([]),[],'Kalman object for state initialization'), ...
        arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level'), ...
        arg({'downsampleFactor','DownsampleFactor'},[],[],'Downsample Factor. Use when memory resources are limited. For k=downsampleFactor, starting from sample t=max(2,k), store only every k Kalman variables (var coefficients (states), etc).'), ...
        arg({'constraints','ConstraintsObject','Constraints'},struct([]),[],'Model constraints. This is a structure with fields .D and .d containing constraints of the form Dx = d (see [1])') ...
        );
    
arg_toworkspace(g);

% check inputs
if isempty(constraints)
    doConstraints = 0;
else
    doConstraints = 1;
    Constr_D = constraints.D;
    Constr_d = constraints.d;
end


[npnts, nchs, ntr] = size(y);
L = nchs*nchs*p;

if npnts<(p+1),
    fprintf(2,'Not enough observed data supplied for given model order\n');
    return
end

% ypred = zeros(size(y));	%prediction of y

% size of downsampled storage
dslen = ceil((npnts-downsampleFactor+(downsampleFactor>1))/downsampleFactor);

VAR=zeros(L,dslen,ntr);

if nargout>1
    PE=zeros(nchs,nchs,dslen,ntr);
end
if nargout>4
    Kout = zeros(nchs,nchs,dslen,ntr);
end;

if verb==2
    h=waitbar(0,sprintf('fitting VAR[%d] model [mode=%s] ...', ...
        p, 'Kalman'));
end

%Kalman Filter initialization (Kp (K predicted or a-priori) equals K(n+1,n) )
F       = eye(L);          % observation matrix
G       = zeros(L,nchs);   % Kalman Gain
x       = zeros(L,1);      % state vector
Kp      = eye(L);
Qstate  = eye(L);          % state noise covariance matrix
Qobs    = eye(nchs);       % measurement noise covariance matrix
ypred   = zeros(size(y));  % prediction of y

%     Kalman=struct('F',eye(L),'H',zeros(nchs,L),'G',zeros(L,nchs),'x',zeros(L,1),'Kp',eye(L),'Qstate',eye(L)*UC,'Qobs',eye(nchs),'ypred',zeros(size(y)));
if ~isempty(Kalman)
    if isfield(Kalman,'ypred'),  ypred = Kalman.ypred;  end
    if isfield(Kalman,'F'),         F  = Kalman.F;      end
    if isfield(Kalman,'Qstate'),Qstate = Kalman.Qstate; end
    if isfield(Kalman,'Kp'),        Kp = Kalman.Kp;     end
    if isfield(Kalman,'Qobs'),    Qobs = Kalman.Qobs;   end
    if isfield(Kalman,'x'),         x  = Kalman.x;      end
    if isfield(Kalman,'H'),         H  = Kalman.H;      end
    if isfield(Kalman,'G'),         G  = Kalman.G;      end
end

upd = eye(L)/L*UC;		%diagonal matrix containing UC

if(mode==3)
    Block=kron(eye(nchs),ones(nchs*p));
elseif(mode==4)
    index=[];
    Block1=[];
    Block0=[];
    for i=1:nchs,
        index=[index ((i-1)*nchs*p+i:nchs:i*nchs*p)];
        mone=eye(nchs);
        mone(i,i)=0;
        mzero=eye(nchs)-mone;
        Block1=blkdiag(Block1,kron(eye(p),mone));
        Block0=blkdiag(Block0,kron(eye(p),mzero));
    end;
elseif mode==5
    Qstate = upd;  % a4 of thesis
elseif mode==6
    Qstate = eye(L)*UC^2;
elseif mode==7
    Qstate = eye(L)*UC;
elseif mode==8
    Qstate = zeros(L);  % RLS algorithm (no process noise)
end


for tr=1:ntr
    
    % TODO: add option to re-initialize the state variables here
    
    curval = 1;
    
    for n = 2:npnts
        
        if verb==2 && ~mod(n,100)
            waitbar(n/npnts,h,...
                {sprintf('Trial (%d/%d)',tr,ntr), ...
                 sprintf('fitting VAR[%d] model [mode=%s] (%d/%d) ...',p,'Kalman',n,npnts)});
        end
                
        if(n<=p)
            Yr=[y(n-1:-1:1,:,tr)' zeros(nchs,p-n+1)];	%vector of past observations
            Yr=Yr(:)';
        else
            Yr=y(n-1:-1:n-p,:,tr)';                     %vector of past observations
            Yr=Yr(:)';
        end
        
        %Update of measurement matrix
        H = blkdiageye(Yr,nchs);
        %     H=kron(eye(nchs),Yr);
        
        %calculate a-priori prediction error (1-step prediction error)
        ypred(n,:,tr)=(H*x)';
        err=y(n,:,tr)-ypred(n,:,tr);
        
        if ~any(isnan(err(:))),
            %update of Qobs (measurement noise covariance matrix, V)) using the prediction error of the previous step
            Qobs=(1-UC)*Qobs+UC*(err'*err);
            
            
            KpH=Kp*H';
            HKp=H*Kp;
            
            %Kalman gain
            G=KpH/(H*KpH+Qobs);
            
            %calculation of the a-posteriori state error covariance matrix
            %K=Kp-G*KpH'; Althouh PK is supposed to be symmetric, this operation makes the filter unstable
            K=Kp-G*HKp;
            
            %mode==0 no update of Qstate (process noise covariance matrix, W)
            %update of Qstate using the predicted state error cov matrix
            if (mode==1)
                Qstate=diag(diag(K)).*UC;
            elseif(mode==2)  % similar to mode a2 from thesis
                Qstate=upd*trace(K);
            elseif(mode==3)
                Qstate=diag(sum((Block*diag(diag(K))),2)')/(p*nchs)*UC;
            elseif(mode==4)
                avg=trace(K(index,index))/(p*nchs)*UC;
                Qstate=Block1*UC+Block0*avg;
            end
            
            %a-priori state error covariance matrix for the next time step
            Kp=K+Qstate;
            
            %current estimation of state x
            x=x+G*(err)';
            
            if doConstraints
                % project the solution onto the constraint surface
                KD = K*Constr_D';
                x = x - (KD/(Constr_D*KD))*(Constr_D*x - Constr_d);
            end
            
        end
        
        if ~mod(n,downsampleFactor)
            % store the current state
            VAR(:,curval,tr) = x;
            
            if nargout>1
                PE(:,:,curval,tr)=Qobs;
            end
            
            if nargout>4
                Kout(:,:,curval,tr) = K;
            end
            
            curval = curval + 1;
        end
    end
    
end

% return residuals
if nargout > 2
    residuals = y - ypred;
end

% return Kalman structure
if nargout > 3
    Kalman.ypred = ypred;
    Kalman.F  = F;
    Kalman.Qstate = Qstate;
    Kalman.Kp = Kp;
    Kalman.Qobs = Qobs;
    Kalman.x = x;
    Kalman.H = H;
    Kalman.G = G;
end

% make sure coefficients are in the
% right order [nchs nchs*p trials]
VAR = permute(VAR,[2 1 3]);

if verb==2, close(h); end
