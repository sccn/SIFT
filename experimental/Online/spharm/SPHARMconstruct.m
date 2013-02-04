function Y = SPHARMconstruct(k,unitSphereFile)
% 
%  SPHARMconstruct=SPHARMconstruct(directory,k)
%
%  directory: directory where SPHARM basis should be saved. For instance,
%                 directory='d:\harmonics\basis\'
%                 WARNINING: The hard drive must have at least 2.4GB of space.
%  k           : degree of spherical harmonics. The recommend value of k is less than 86 otherwise 
%                 you will encounter the numerical singularity caused by MATLAB. 
%
%  Consruct spherical harmonic functions upto the degree k. 
%  Requires two external files Y_l.m and unitsphere.mat
%
% (C) Moo K. Chung, 2006-2008
%       Department of Biostatistics and Medical Informatics
%       University of Wisconsin-Maison
%  
% email://mkchung@wisc.edu
% http://www.stat.wisc.edu/softwares/weighted-SPHARM/weighted-SPHARM.html
%
% If you use this code, please reference the following paper. 
% You need to read the paper to modify the code.
%
% Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007. 
% Weighted Fourier series representation and its application to quantifying 
% the amount of gray matter. IEEE Transactions on Medical Imaging, 26:566-581.
%
% Update history: 
% Created Sept 19 2006; Modified July 5, 2007 
%-------------------------------------------------------------------------------
directory = fileparts(which('SPHARMconstruct'));
if ~exist([directory filesep 'basis'],'dir'), mkdir([directory filesep 'basis']);end
directory = [directory filesep 'basis' filesep];
if nargin < 2
    load unitsphere.mat;
else
    [coord,tri] = freesurfer_read_surf(unitSphereFile);
    coord = coord';
    [theta,varphi,nbr] = cart2sph(coord(1,:),coord(2,:),coord(3,:));
end
% load unitsphere_reduced
Y = [];
%tic
hwait = waitbar(0,'Constructing SH basis...');
for l=0:k
    %elapsedtime(l+1)=toc; these lines for measuring running time
    %[l elapsedtime(l+1)]
    f_name = [directory int2str(l) '.mat'];
    if ~exist(f_name,'file')
        temp   = Y_l(l,theta,varphi);
    else
        load(f_name);
    end
    if nargout
        Y = [Y; temp];
    else
        save(f_name,'temp')
    end
    waitbar(l/k,hwait);
end;
close(hwait);
