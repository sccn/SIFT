% SPHARMconstruct(85);
% figure;trisurf(tri,coord(1,:),coord(2,:),coord(3,:))
clear all;clc;close all;

surfacefile = '/media/DATOS/my_work/PhilTrans/smoothed_surfaces/corrected_MC0000005_source_BEM_mesh.mat';
load(surfacefile);
patch(SurfData,'facecolor','g')

surfacefile = '/home/ale/MNI/Linux-i686/CIVET/models/surf_reg_model_right.obj';
[tri,coord,nbr,normal]=mni_getmesh('outersurface.obj');
fv.faces = tri;
fv.vertices = coord';
patch(fv,'facecolor','g')

L=40;
sigma=0;%0.0001;
[coord_smooth,fourier_coeff]= SPHARMsmooth(fv.vertices',L,sigma);
% [coord_smooth2,fourier_coeff]= SPHARMsmooth(fv.vertices(end/2+1:end,:)',L,sigma);
fv2 = fv;
fv2.vertices = coord_smooth';
figure;patch(fv2,'facecolor','g')


[fv.vertices,fv.faces] = getFV(4);
figure;patch(fv)