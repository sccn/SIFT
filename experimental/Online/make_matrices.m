clear all
clc
close all
folder = '/home/alejandro/testing_connectome';

[Surface.lh.vertices,Surface.lh.faces] = freesurfer_read_surf([folder filesep 'lh.pial']);
[Surface.rh.vertices,Surface.rh.faces] = freesurfer_read_surf([folder filesep 'rh.pial']);

[~, Surface.lh.label, Surface.lh.colortable] = read_annotation([folder filesep 'lh.aparc.annot']);
[~, Surface.rh.label, Surface.rh.colortable] = read_annotation([folder filesep 'rh.aparc.annot']);

for it=1:size(Surface.lh.colortable.table,1)
    Surface.lh.colortable.struct_names{it} = ['left ' lower(Surface.lh.colortable.struct_names{it})];
    Surface.rh.colortable.struct_names{it} = ['right ' lower(Surface.rh.colortable.struct_names{it})];
end

nodes = load([folder filesep 'nodes.txt']);
start_pos = load([folder filesep 'start_pos'])';
end_pos = load([folder filesep 'end_pos'])';
fid = fopen([folder filesep 'regions.txt'],'r');
it =1;
count = 1;
deep_region = [];
while ~feof(fid)
    tmp = lower(fgets(fid));
    if str2double(sprintf('%i',tmp(end))) == 10, tmp(end) = [];end
    Regions{count,1} = tmp; %#ok
    if ~isempty(strfind(tmp,'ctx'))
        ind = find(tmp == '-');
        if length(ind) >= 2
            if strcmp(tmp(ind(1)+1:ind(2)-1),'rh')
                regions{it,1} = ['right ' tmp(ind(2)+1:end)]; %#ok
            elseif strcmp(tmp(ind(1)+1:ind(2)-1),'lh')
                regions{it,1} = ['left ' tmp(ind(2)+1:end)];  %#ok
            else
                tmp(ind(1)) = ' ' ;
                regions{it,1} = tmp;    %#ok
            end
        elseif length(ind) == 1
            tmp(ind(1)) = ' ' ;
            regions{it,1} = tmp;    %#ok
        else
            regions{it,1} = tmp;                              %#ok
        end
        regions{it,1}(regions{it,1}=='-') = []; %#ok
        it = it+1;
    else
        deep_region(end+1) = count;%#ok
    end
    count = count+1;
end
fclose(fid);

A = cat(1,Surface.lh.colortable.struct_names,Surface.rh.colortable.struct_names);
rm_structure = setdiff(A,regions);
Surface.lh = removeStructure(Surface.lh,rm_structure(1));
Surface.rh = removeStructure(Surface.rh,rm_structure(3));

[SurfData.vertices,SurfData.faces] = mergemesh(Surface.lh.vertices,Surface.lh.faces,Surface.rh.vertices,Surface.rh.faces);
SurfData.faces(:,4) = [];
[SurfData.vertices,SurfData.faces] = meshcheckrepair(SurfData.vertices,SurfData.faces,'duplicated');
[SurfData.vertices,SurfData.faces] = meshcheckrepair(SurfData.vertices,SurfData.faces,'isolated');
[SurfData.vertices,SurfData.faces] = meshcheckrepair(SurfData.vertices,SurfData.faces,'deep');

figure;plot_over_surf(SurfData);camlight('headlight')
SurfData.label = [Surface.lh.label; max(Surface.lh.label) + Surface.rh.label];
SurfData.colortable = Surface.lh.colortable;
SurfData.colortable.numEntries = Surface.lh.colortable.numEntries + Surface.rh.colortable.numEntries;
SurfData.colortable.struct_names = cat(1,Surface.lh.colortable.struct_names,Surface.rh.colortable.struct_names);
SurfData.colortable.table = [Surface.lh.colortable.table; Surface.rh.colortable.table];
SurfData.colortable.table(:,5) = [Surface.lh.colortable.table(:,5); max(Surface.lh.label)+Surface.rh.colortable.table(:,5)];

figure;
plotSurf(SurfData)
camlight('headlight')

nodes(deep_region,:) = [];

[~,loc] = ismember(start_pos,nodes,'rows');
loc = ~loc;
start_pos(loc,:) = [];
end_pos(loc,:) = [];

[~,loc] = ismember(end_pos,nodes,'rows');
loc = ~loc;
start_pos(loc,:) = [];
end_pos(loc,:) = [];

loc = all(start_pos == end_pos,2);
start_pos(loc,:) = [];
end_pos(loc,:) = [];

[~,I] = ismember(start_pos,nodes,'rows');
[~,J] = ismember(end_pos,nodes,'rows');
C = zeros(length(regions));
for it=1:length(I), C(I(it),J(it)) = 1;C(J(it),I(it)) = 1;end
plotConnectivityMatrix(C,regions);

anatomicalModel = plotModel(SurfData,C,regions);

% resampling source space
decimationRate = 0.01; % reducing the mesh to its 1%
[rSurfData.vertices,rSurfData.faces]=meshresample(SurfData.vertices,SurfData.faces,decimationRate);

figure; patch('vertices',SurfData.vertices,'faces',SurfData.faces,'FaceColor','g','facelighting','gouraud','LineStyle','-','facealpha',0.75);
hold on;patch('vertices',rSurfData.vertices,'faces',rSurfData.faces,'FaceColor',[0.9922 0.9176 0.7961],'facelighting','gouraud','LineStyle','-.');

F = TriScatteredInterp(SurfData.vertices(:,1),SurfData.vertices(:,2),SurfData.vertices(:,3),SurfData.label,'nearest');
rSurfData.label = F(rSurfData.vertices(:,1),rSurfData.vertices(:,2),rSurfData.vertices(:,3));
rSurfData.colortable = SurfData.colortable;


rAnatomicalModel = plotModel(rSurfData,C,regions);
save('dataset.mat','anatomicalModel','rAnatomicalModel');