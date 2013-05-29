function disp=SPHARMcorrespondence(coeff1,coeff2,sigma)
%---------------------------------------------------------------------------------------
%disp=SPHARMcorrespondence(coeff1,coeff2)
%
%coeff1    :     SPHARM coefficients of surface 1
%coeff2    :     SPHARM coefficients of surface 2
%                coeff1 (and also coeff 2) is a structured array containg coeff1.x, coeff1.y, coeff1.z
%                For degree k representation, coeff.x is the SPHARM coefficients of x-cooridinate given 
%                as (k+1)*(2k+1) matrix. At degree k, there are 2k+1 order harmonics.
%                Ex. coeff.x(2+1,1:(2*2+1)) is the 2nd degree SPHARM coefficients 
%                of all orders -2, -1, 0, 1, 2.  
%sigma      :    bandwidth. it controls the amount of smoothing. Note that sigma is not FWHM. 
%                See the paper below for detail. 
%                
%disp      :     matrix of size 3 * n_vertices respresenting the diplacement vector 
%                field from surface 1 to surface 2. disp(2,34) is the y-coordinate of 
%                the vector for 34-th mesh vertex. 
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
% Update history: Created July 6, 2007
%----------------------------------------------------------------------------------------


% The algorithm assumes there are n_vertex = 40962 on unit sphere which servers as a grid.
directory = 'd:\harmonics\basis\';   %diretory where SPHARM bases are stored. 
n_vertex = 40962   
% Initialization for displacement vector
xcoord=zeros(n_vertex,1);
ycoord=zeros(n_vertex,1);
zcoord=zeros(n_vertex,1);

% Determine degree k
dim=size(coeff1.x);
k=dim-1; 

% 0-th degree (l=0). This initial condition doesn't fit into the iteration below.
    f_name=strcat(directory,int2str(0));
    load(f_name);
    Y=temp';
    clear temp;
    betal= coeff2.x(1,1)-coeff1.x(1,1);
    xcoord =Y*betal';
    betal= coeff2.y(1,1)-coeff1.y(1,1);
    ycoord=Y*betal';
    betal= coeff2.z(1,1)-coeff1.z(1,1);
    zcoord=Y*betal';
        
for l=1:k
    l   
      
    f_name=strcat(directory,int2str(l));
    load(f_name);
    Y=temp';
    % real(Y_l) gives negative harmonics imag(Y_l) gives positive
    % harmonics 
    clear temp;
    Y=[real(Y)   imag(Y(:,2:(l+1)))];
    
    betal= coeff2.x(l+1,1:2*l+1)-coeff1.x(l+1,1:2*l+1);
    xcoord =xcoord + exp(-l*(l+1)*sigma)*Y*betal';
    
    betal= coeff2.y(l+1,1:2*l+1)-coeff1.y(l+1,1:2*l+1);
    ycoord=ycoord + exp(-l*(l+1)*sigma)*Y*betal';
    
    betal= coeff2.z(l+1,1:2*l+1)-coeff1.z(l+1,1:2*l+1);
    zcoord=zcoord+ exp(-l*(l+1)*sigma)*Y*betal';
    
end;

disp=[xcoord, ycoord, zcoord]';
return;
