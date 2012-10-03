load /Users/timmullen/Documents/WORK/SIFT/Code/SIFT_bitbucket/vis/@resources/LONImesh_fullres.mat

figure; 

hlp_plotAtlas(lonimesh);
box off
axis off
axis vis3d
camlight left
camlight right
lighting phong
material metal

surfh = findobj(get(findobj(get(gcf,'children'),'type','axes'),'children'),'type','patch');
nvertices = size(lonimesh.vertices,1);

%% brainbow mode
set(surfh,'FaceVertexCData',rainbow(nvertices))

%% Gazzaniga mode
cmap = jet(nvertices);
% for i=1:3
%     cmapg(:,i) = interleave(cmap(1:end/2,i), cmap(end/2+1:end,i));
% end

cmapg = cmap;

set(surfh,'FaceVertexCData',cmapg);

set(surfh,'FaceAlpha','flat');

set(surfh,'FaceVertexAlphaData',[zeros(nvertices/2,1); ones(nvertices/2,1)]);   %interleave()'
