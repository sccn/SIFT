function [coord_smooth, coeff]=SPHARMsmooth(coord,L,sigma,n_vertex)
%---------------------------------------------------------------------------------------
%[coord_smooth, coeff]=SPHARMsmooth(coord,L,sigma)
%
%
% coord            : Either 2D or 3D matrix of coordinates of cortical surface. 
%                    For 3D, the dimension of matrix is n_subject * 3 cooridnates * n_vertex.
%                    For 2D, the dimension of matrix is 3 coordinates * n_vertex.  
%                    For 3D matrix, coord(1,2,:) is the all the y-coordinates of the first subject. 
%                    For MNI mesh, we have n_vertex=40962.
%
%    L               : The maximal degree of weighted-SPHARM representation.
%                       Read the paper below to find it optimally.
%
%  sigma          : bandwith of weighted-SPHARM representation
%                       It is the bandwidth of heat kernel on unit sphere.
%                       When sigma=0, it is the traditional SPHARM representation. 
%                       range beween 0.0001 and 0.01 will be sufficient for cortical surfaces.
%
% coord_smooth: The weighted-SPHARM result. The dimension is identical to coord.
%
% coeff   : The estimated SPHARM coefficients (Fourier coefficients) given as a structured array
%                        containg coeff.x, coeff.y, coeff.z
%                        coeff.x is the SPHARM coefficients of x-cooridinate given as (L+1)*(2L+1)*n_subject 
%                        coeff.x(3,:,10) is the 2nd degree SPHARM coefficients of all order for the 10th subject.  
%
%
%
% (C) Moo K. Chung, 2006, 2007
%     Department of Biostatistics and Medical Informatics
%     Waisman Laboratory for Brain Imaging
%     University of Wisconsin-Maison
%  
% email://mchung@stat.wisc.edu
% http://www.stat.wisc.edu/softwares/weighted-SPHARM/weighted-SPHARM.html
%
% If you use this code, please reference the following paper. 
% You need to read the paper to understand the notations and the algorithm.
%
% Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007. 
% Weighted Fourier series representation and its application to quantifying 
% the amount of gray matter. IEEE Transactions on Medical Imaging, 26:566-581.
%
%
% Update history Sept 19 2006; July 5, 2007
%----------------------------------------------------------------------------------------

if nargin<4
    n_vertex = 40962;   % the number of vertices in a mesh. MNI mesh has the fixed number of vertices for all subjects.
end
directory = [fileparts(which('SPHARMconstruct')) filesep 'basis' filesep];   %diretory where SPHARM bases are stored. 
                                                       

% determining the number of subjects
if ndims(coord)==3
    n_subject=size(coord,1);
else
    n_subject=1;
end;

% initialization
if n_subject==1
    xestimate=zeros(n_vertex,1);
    yestimate=zeros(n_vertex,1);
    zestimate=zeros(n_vertex,1);
else
    xestimate=zeros(n_vertex,n_subject);
    yestimate=zeros(n_vertex,n_subject);
    zestimate=zeros(n_vertex,n_subject);
end;

betax=zeros(L+1,2*L+1,n_subject);
betay=zeros(L+1,2*L+1,n_subject);
betaz=zeros(L+1,2*L+1,n_subject);

% 0-degree SPHARM coefficients. Step 2 in the iterative resiual fitting (IRF) algorithm. See the paper.
for i=1:n_subject
    
    if n_subject==1
        x=coord(1,:)';
        y=coord(2,:)';
        z=coord(3,:)';
    else
        x=squeeze(coord(i,1,:));
        y=squeeze(coord(i,2,:));
        z=squeeze(coord(i,3,:));
    end;
  
    f_name=[directory int2str(0) '.mat'];
    load(f_name);
    Y=temp';
    Ycommon=(Y'*Y)\Y';
    
    betal=Ycommon*x;
    betax(1,1,i)=betal';
    xsmooth=Y*betal;
    xestimate(:,i)=xsmooth;
    
    betal=Ycommon*y;
    betay(1,1,i)=betal';
    ysmooth=Y*betal;
    yestimate(:,i)=ysmooth;
    
    betal=Ycommon*z;
    betaz(1,1,i)=betal';
    zsmooth=Y*betal;
    zestimate(:,i)=zsmooth;
end;


hwait = waitbar(0,'Smoothing...');
for l=1:L
    % degree l SPHARM bases are recycled for all subject reducing computational time.
    
    for i=1:n_subject
        if n_subject==1
            x=coord(1,:)';
            y=coord(2,:)';
            z=coord(3,:)';
        else
            x=squeeze(coord(i,1,:));
            y=squeeze(coord(i,2,:));
            z=squeeze(coord(i,3,:));
        end;
        
        % Step 4: residual. See the paper for detail 
        x_j = x-xestimate(:,i);
        y_j = y-yestimate(:,i);
        z_j = z-zestimate(:,i);
   
        f_name=strcat(directory,int2str(l));
        load(f_name);
        Y=temp';
        % real(Y_l) gives negative harmonics imag(Y_l) gives positive harmonics 
        clear temp;
        Y=[real(Y)   imag(Y(:,2:(l+1)))];
        Ycommon=(Y'*Y)\Y';
        
        % Step 5: refitting the residual. See the paper for detail
        betal=Ycommon*x_j;
        betax(l+1,1:2*l+1,i)=betal';
        xsmooth=Y*betal;
        xestimate(:,i)=xestimate(:,i)+ exp(-l*(l+1)*sigma)*xsmooth;
        
        betal=Ycommon*y_j;
        betay(l+1,1:2*l+1,i)=betal';
        ysmooth=Y*betal;
        yestimate(:,i)=yestimate(:,i)+ exp(-l*(l+1)*sigma)*ysmooth;
        
        betal=Ycommon*z_j;
        betaz(l+1,1:2*l+1,i)=betal';
        zsmooth=Y*betal;
        zestimate(:,i)=zestimate(:,i)+ exp(-l*(l+1)*sigma)*zsmooth;
    end;
    waitbar(l/L,hwait);
end;
close(hwait);

%output the results in a proper shape
coord_smooth=[xestimate yestimate zestimate]';
coord_smooth=squeeze(reshape(coord_smooth,n_subject,3,n_vertex));

coeff.x=betax;
coeff.y=betay;
coeff.z=betaz;
return;