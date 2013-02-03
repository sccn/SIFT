function fourier_vector=SPHARMvectorize(fourier_coeff, k)
%
% It vectorizes the fourier coefficient matrix and removes the padded zeros up to degree k.
%---------------------------------------------------------------------------------------
%
% (C) Moo K. Chung, 2008
%     Department of Biostatistics and Medical Informatics
%     Waisman Laboratory for Brain Imaging
%     University of Wisconsin-Maison
%  
% email://mkchung@wisc.edu
% http://www.stat.wisc.edu/softwares/weighted-SPHARM/weighted-SPHARM.html
%
% If you use this code, please reference the following paper. 
% You need to read the paper to understand the notations and the algorithm.
%
% Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007. 
% Weighted Fourier series representation and its application to quantifying 
% the amount of gray matter. IEEE Transactions on Medical Imaging, 26:566-581.
%
% Update history: Created January 6, 2008
%----------------------------------------------------------------------------------------

s=size(fourier_coeff.x);
subject=s(3);
degree=k+1;
for i=1:subject
    x=fourier_coeff.x(:,:,i);
    y=fourier_coeff.y(:,:,i);
    z=fourier_coeff.z(:,:,i);
    vec_x=x(1,1);
    vec_y=y(1,1);
    vec_z=z(1,1);
    for j=2:degree
        vec_x=[vec_x x(j,1:(2*j-1))];
        vec_y=[vec_y y(j,1:(2*j-1))];
        vec_z=[vec_z z(j,1:(2*j-1))];
    end

    fourier_vector.x(:,i)=vec_x;
    fourier_vector.y(:,i)=vec_y;
    fourier_vector.z(:,i)=vec_z;
end;


