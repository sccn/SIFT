load LONImesh_fullres_new.mat
metadata = importdata('LONI_anatomy.txt',' ');

lblidx = cell2mat(cellfun(@str2num,metadata.textdata(:,1),'UniformOutput',false));

lonimesh       = Surf;
lonimesh.label = filtLabels;

n = length(lblidx);     % number of unique labels (anatomical regions)
cmap = distinguishable_colors(n,[1 1 1; 0 0 0]);

lonimesh.colortable.numEntries = n;
lonimesh.colortable.orig_tab = '';
lonimesh.colortable.struct_names = metadata.textdata(:,2);
lonimesh.colortable.table = [cmap, zeros(n,1), lblidx];

       

