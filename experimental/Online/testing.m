%SPHARMconstruct(85);
clear all;clc;close all;
load /home/alejandro/testing_connectome/dataset.mat
[fv.faces,fv.vertices]=mni_getmesh('outersurface.obj');
fv.vertices = fv.vertices';
figure;
colordef black
patch(fv,'linestyle','none','facealpha',0.075)
hold on
line([start_pos(:,1) end_pos(:,1)],[start_pos(:,2) end_pos(:,2)],[start_pos(:,3) end_pos(:,3)])
scatter3(nodes(:,1),nodes(:,2),nodes(:,3),'filled');
for it=1:length(regions), text('Position',nodes(it,:),'String',regions{it});end




%%

L=40;
sigma=0;%0.0001;
[Surface.vertices,fourier_coeff]= SPHARMsmooth(Surface.vertices',L,sigma,size(Surface.vertices,1));
% [coord_smooth2,fourier_coeff]= SPHARMsmooth(fv.vertices(end/2+1:end,:)',L,sigma);
Surface.vertices = Surface.vertices';


% 
% [Surface.vertices,Surface.faces] = getFV(4);
% figure;patch(fv)



%%

SPHARMconstruct(85);
