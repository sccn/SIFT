function [tri,coord,nbr,normal]=mni_getmesh(inputmesh)
% [tri,coord,nbr,normal]=mni_getmesh(inputmesh)
%
% tri    : list of triangle elements
% coord  : node coordinatates 
% nbr    : 1st neihbor nodes list 
% normal : normal vector list
%
% EXAMPLE: [tri,coord,nbr,normal]=mni_getmesh(inputmesh)
%
%(C) Moo K. Chung
% Update history Feb 5. 2004; July 20, 2004; Oct 28, 2004; Sept 7, 2007
% mkchung@wisc.edu
% http://www.stat.wisc.edu/softwares/hk/hk.html

% This is the montreal neurological institute (MNI) specific ASCII triangular mesh data structure.
% For FreeSurfer software, a slightly different data input coding is needed. It will be provided upon request.
%
% Mesh topology:
% There are total V + F rows.
% V = number of vertices
% F = number of faces
% It is assumed that the input mesh is topologically equvalent to a sphere, i.e.
% F = 2V - 4: Euler characteristic equation.
% Running time = 30 seconds in Pentinum-M 1.5Ghz.


fid=fopen(inputmesh);
frewind(fid);
fscanf(fid,'%c',1)
fscanf(fid,'%f',5)
n_points=fscanf(fid,'%i',1)
coord=fscanf(fid,'%f',[3,n_points]);
normal=fscanf(fid,'%f',[3,n_points]);
n_tri=fscanf(fid,'%i',1)
fscanf(fid,'%i',5+n_tri);
tri=fscanf(fid,'%i',[3,n_tri])'+1;
fclose(fid);



% compute the maximum degree of node
degree=zeros(n_points,1);
for j=1:n_tri
degree(tri(j,:))=degree(tri(j,:))+1;
end
max_degree=max(degree);


% find out the neighboring nodes
nbr=zeros(n_points,max_degree);
for i_tri=1:n_tri
    for j=1:3
        cur_point = tri(i_tri,j);
        for k=1:3
            if (j ~= k)
                nbr_point= tri(i_tri,k);
                if find(nbr(cur_point,:)==nbr_point)
                    ;
                else
                    n_nbr = min(find(nbr(cur_point,:) == 0));
                    nbr(cur_point,n_nbr) = nbr_point;
                end;
            end;
        end;
    end;
end;

return;


